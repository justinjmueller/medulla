"""
Abstract base class and concrete implementations for dataset definitions
used in the SPINE production pipeline.
"""
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional
import logging
import samweb_client as sam
from glob import glob
from tqdm import tqdm
from pathlib import Path
import pandas as pd
import uproot
import h5py
import sqlite3

from library import configure_samweb, command

logger = logging.getLogger(__name__)


class Dataset(ABC):
    """
    Abstract base class for CAF, HDF5, and Match datasets used in the
    SPINE production pipeline.
    """
    def __init__(
        self,
        name: str,
        project: str,
        strategies: list[str],
        path: dict[str, str],
        register: Path,
    ):
        self.name = name
        self.project = project
        self.strategies = strategies
        self.path = path
        self.register = register
        self.data: Optional[pd.DataFrame] = None

    @abstractmethod
    def process(self) -> None:
        pass

    def process_strategy_load(self) -> bool:
        """Load dataset metadata from a pre-existing CSV in the register."""
        csv_file = self.register / (self.name + '.csv')
        if csv_file.exists():
            self.data = pd.read_csv(csv_file)
            return True
        return False

    def write(
        self,
        output_path: Path
    ) -> None:
        """
        Write dataset metadata to a CSV file. Skips the write if the
        existing file already matches the current data.
        """
        if self.data is None:
            logger.warning(
                f"No data available for dataset {self.name}; "
                f"skipping write."
            )
            return

        if output_path.exists():
            existing_data = pd.read_csv(output_path)
            if existing_data.equals(self.data):
                logger.info(
                    f"Dataset {self.name} at {output_path} is up to "
                    f"date; skipping write."
                )
                return

        self.data.to_csv(output_path, index=False)
        logger.info(f"Wrote dataset {self.name} to {output_path}.")


class CAFDataset(Dataset):
    """
    Concrete implementation of a CAF dataset. Supports three file
    discovery strategies: load (pre-existing CSV), sam (SAMWeb query),
    and manual (filesystem glob + uproot metadata extraction).
    """

    def __init__(
        self,
        *args,
        name: str,
        project: str,
        strategies: list[str],
        path: dict[str, str],
        register: Path,
        **kwargs,
    ):
        super().__init__(
            name=name,
            project=project,
            strategies=strategies,
            path=path,
            register=register,
        )

    def process(
        self,
        experiment: str
    ) -> None:
        """
        Iterate through strategies in order until one succeeds.
        """
        success = False
        for strategy in self.strategies:
            if strategy == 'load':
                success = self.process_strategy_load()
            elif strategy == 'sam':
                samweb = configure_samweb(exp=experiment)
                success = self.process_strategy_sam(samweb)
            elif strategy == 'manual':
                success = self.process_strategy_manual()
            if success:
                break

        logger.info(
            f"Completed processing for CAF dataset {self.name} "
            f"({len(self.data) if self.data is not None else 0} rows)."
        )

    @staticmethod
    def _parse_event_list(raw: str) -> list:
        """
        Parse a SAMWeb event-number list string into a sorted list of ints.

        Data CAF files store events as an underscore-delimited string with a
        leading underscore, e.g. ``_49208_49340_49780_``.  Returns an empty
        list if the string is absent or malformed.
        """
        if not raw:
            return []
        return [int(e) for e in raw.strip('_').split('_') if e]

    def process_strategy_sam(
        self,
        samweb: sam.SAMWebClient
    ) -> bool:
        """
        Query SAMWeb for dataset files and extract run/event metadata.

        Handles both simulation (file_type='mc') and data (file_type='data').
        Simulation files expose top-level ``first_event`` and ``event_count``
        fields.  Data files carry the equivalent information in the custom
        parameters ``sbnd.event_number_list`` (first event) and
        ``sbn_dm.event_count`` (total events).
        """
        meta = {
            'file': [],
            'run': [],
            'subrun': [],
            'first_event': [],
            'total_events': []
        }
        try:
            files = samweb.listFiles(defname=self.path['sam'])
            if not files:
                return False

            progress = tqdm(
                range(0, len(files), 1000),
                desc="Looking up metadata",
                unit="batch"
            )
            for i in progress:
                chunk = files[i:i+1000]
                res = samweb.getMultipleMetadata(chunk, locations=True)
                for x in res:
                    loc = x['locations'][0]['location'].split(':', 1)[-1]
                    meta['file'].append(loc + '/' + x['file_name'])
                    meta['run'].append(x['runs'][0][0])
                    meta['subrun'].append(x['runs'][0][1])

                    if x.get('file_type') == 'data':
                        events = self._parse_event_list(
                            x.get('sbnd.event_number_list', '')
                        )
                        meta['first_event'].append(events[0] if events else 0)
                        meta['total_events'].append(
                            int(x.get('sbn_dm.event_count', len(events)))
                        )
                    else:
                        meta['first_event'].append(x['first_event'])
                        meta['total_events'].append(x.get('event_count', 0))

            self.data = pd.DataFrame(meta).sort_values(
                by=['run', 'first_event']
            )
            return True

        except Exception as e:
            logger.exception(
                f"Error processing SAM dataset {self.name}: {e}"
            )
            return False

    def process_strategy_manual(self) -> bool:
        """
        Discover CAF files via filesystem glob and extract metadata with
        uproot.
        """
        keys = [
            'rec/hdr/hdr.run',
            'rec/hdr/hdr.subrun',
            'rec/hdr/hdr.evt',
        ]
        meta = {
            'file': [],
            'run': [],
            'subrun': [],
            'first_event': [],
            'total_events': []
        }
        try:
            pattern = self.path['pattern']
            files = glob(pattern)
            files = [f for f in files if 'flat.' not in f]
            if not files:
                return False

            progress = tqdm(files, desc="Extracting metadata", unit="file")
            for f in progress:
                with uproot.open(f) as uf:
                    if 'recTree' not in uf:
                        meta['file'].append(f)
                        meta['run'].append(None)
                        meta['subrun'].append(None)
                        meta['first_event'].append(None)
                        meta['total_events'].append(0)
                        continue

                    info = uf['recTree'].arrays(keys, library='pd')
                    meta['file'].append(f)
                    meta['run'].append(info['rec/hdr/hdr.run'].iloc[0])
                    meta['subrun'].append(info['rec/hdr/hdr.subrun'].iloc[0])
                    meta['first_event'].append(info['rec/hdr/hdr.evt'].iloc[0])
                    meta['total_events'].append(len(info))

            self.data = pd.DataFrame(meta).sort_values(
                by=['run', 'first_event']
            )
            return True

        except Exception as e:
            logger.exception(
                f"Error processing manual CAF dataset {self.name}: {e}"
            )
            return False


class HDF5Dataset(Dataset):
    """
    Concrete implementation of an HDF5 dataset representing SPINE
    reconstruction outputs. Supports load and manual strategies.
    """

    def __init__(
        self,
        *args,
        name: str,
        project: str,
        strategies: list[str],
        path: dict[str, str],
        register: Path,
        **kwargs,
    ):
        super().__init__(
            name=name,
            project=project,
            strategies=strategies,
            path=path,
            register=register,
        )

    def process(
        self,
        experiment: str
    ) -> None:
        """
        Iterate through strategies in order until one succeeds.
        """
        success = False
        for strategy in self.strategies:
            if strategy == 'load':
                success = self.process_strategy_load()
            elif strategy == 'sam':
                logger.warning("HDF5 SAM strategy not yet implemented.")
            elif strategy == 'manual':
                success = self.process_strategy_manual()
            if success:
                break

        logger.info(
            f"Completed processing for HDF5 dataset {self.name} "
            f"({len(self.data) if self.data is not None else 0} rows)."
        )

    def process_strategy_manual(self, workers: int = 16) -> bool:
        """
        Discover HDF5 files via filesystem glob and extract run_info
        metadata in parallel.

        Reads only the first row of run_info (for run/subrun/first_event)
        and uses the dataset shape attribute for total_events, avoiding
        loading the full dataset into memory.

        Parameters
        ----------
        workers : int
            Number of parallel threads for file I/O. Default is 16,
            which works well for PNFS/dCache.
        """
        pattern = self.path['pattern']
        files = glob(pattern)
        logger.debug(f"HDF5 manual pattern: {pattern}")
        logger.debug(f"HDF5 manual files found: {len(files)}")
        if not files:
            return False

        def _read_meta(f: str) -> tuple:
            """Return (file, run, subrun, first_event, total_events)."""
            with h5py.File(f, 'r') as hf:
                ds = hf['run_info']
                n = ds.shape[0]
                if n == 0:
                    return (f, None, None, None, 0)
                row = ds[0]
                return (f, int(row[0]), int(row[1]), int(row[2]), n)

        rows = []
        n_workers = min(workers, len(files))
        try:
            with ThreadPoolExecutor(max_workers=n_workers) as pool:
                futures = {pool.submit(_read_meta, f): f for f in files}
                progress = tqdm(
                    as_completed(futures),
                    total=len(futures),
                    desc="Extracting metadata",
                    unit="file"
                )
                for future in progress:
                    try:
                        rows.append(future.result())
                    except Exception as e:
                        logger.warning(
                            f"Failed to read {futures[future]}: {e}"
                        )

        except Exception as e:
            logger.exception(
                f"Error processing manual HDF5 dataset {self.name}: {e}"
            )
            return False

        if not rows:
            return False

        self.data = pd.DataFrame(
            rows,
            columns=['file', 'run', 'subrun', 'first_event', 'total_events']
        ).sort_values(by=['run', 'first_event'])
        return True


class MatchDataset(Dataset):
    """
    Concrete implementation of a match between an HDF5 dataset and a
    CAF dataset. A successful match produces the input to the grid-based
    merging step of the SPINE production pipeline.
    """

    def __init__(
        self,
        *args,
        name: str,
        project: str,
        strategies: list[str],
        depends_on: list[str],
        path: dict[str, str],
        register: Path,
        projects_dir: Path,
        **kwargs,
    ):
        super().__init__(
            name=name,
            project=project,
            strategies=strategies,
            path=path,
            register=register,
        )
        self.depends_on = depends_on
        self.projects_dir = projects_dir

    def process_strategy_load(self) -> bool:
        """Load match data from CSV; treats an empty file as a cache miss."""
        if super().process_strategy_load():
            if len(self.data) > 0:
                return True
            self.data = None
        return False

    def process(
        self,
        samples: dict,
        experiment: str,
    ) -> None:
        """
        Iterate through strategies in order until one succeeds.
        `create_project()` is called whenever the project DB does not yet
        exist, regardless of whether data came from the match or load strategy.
        """
        for strategy in self.strategies:
            if strategy == 'load':
                if self.process_strategy_load():
                    break
            elif strategy == 'match':
                if self.process_strategy_match(samples, experiment):
                    break

        logger.info(
            f"Completed processing for Match dataset {self.name} "
            f"({len(self.data) if self.data is not None else 0} rows)."
        )

        project_path = self.projects_dir / ('proj_' + self.name + '.db')
        if not project_path.exists():
            self.create_project()

    def process_strategy_match(
        self,
        samples: dict,
        experiment: str,
    ) -> bool:
        """
        Load the dependent CAF and HDF5 datasets and inner-join them on
        (run, subrun, first_event, total_events).
        """
        if len(self.depends_on) != 2:
            logger.error(
                f"Match dataset {self.name} must depend on exactly two "
                f"datasets."
            )
            return False

        dep_datasets = []
        for dep_name in self.depends_on:
            if dep_name not in samples:
                logger.error(
                    f"Dependent dataset {dep_name} not found for match "
                    f"dataset {self.name}."
                )
                return False
            dep_datasets.append(samples[dep_name])

        caf_dataset = None
        hdf5_dataset = None
        for ds in dep_datasets:
            if ds['type'] == 'caf':
                caf_dataset = ds
            elif ds['type'] == 'spine':
                hdf5_dataset = ds
        if caf_dataset is None or hdf5_dataset is None:
            logger.error(
                f"Match dataset {self.name} must depend on one CAF "
                f"dataset and one HDF5 dataset."
            )
            return False

        caf = CAFDataset(**caf_dataset, register=self.register)
        caf.process(experiment=experiment)

        hdf5 = HDF5Dataset(**hdf5_dataset, register=self.register)
        hdf5.process(experiment=experiment)

        matches = pd.merge(
            caf.data,
            hdf5.data,
            on=['run', 'subrun', 'first_event', 'total_events'],
            suffixes=('_caf', '_hdf5')
        )
        self.data = matches.sort_values(by=['run', 'first_event'])
        return True

    def create_project(self) -> None:
        """
        Write the matched pairs into a SQLite project database. Each row
        represents one grid job (one CAF + one HDF5 file to merge).
        """
        SCHEMA = """
        CREATE TABLE IF NOT EXISTS jobs (
            jobid INTEGER PRIMARY KEY,
            caf_file TEXT NOT NULL,
            hdf5_file TEXT NOT NULL,
            status TEXT NOT NULL DEFAULT 'pending',
            val_path TEXT,
            output_path TEXT
        );
        """

        if self.data is None:
            logger.warning(
                f"No data available for match dataset {self.name}; "
                f"skipping project creation."
            )
            return

        project_path = self.projects_dir / ('proj_' + self.name + '.db')
        project_path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(project_path)
        curs = conn.cursor()

        command(curs, SCHEMA)

        ins = []
        for jobid, (_, row) in enumerate(self.data.iterrows()):
            ins.append((
                jobid,
                row['file_caf'],
                row['file_hdf5'],
                'pending'
            ))
        command(
            curs,
            "INSERT INTO jobs (jobid, caf_file, hdf5_file, status)"
            " VALUES (?, ?, ?, ?)",
            ins
        )
        conn.commit()
        conn.close()

        logger.info(
            f"Wrote project file for match dataset {self.name} to "
            f"{project_path}. Entered {len(ins)} jobs."
        )
