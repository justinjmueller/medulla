/**
 * @file spinevar.h
 * @brief Implementation of functions that create SpillMultiVar objects from
 * constituent functions that implement variables and cuts.
 * @details This file contains the implementation of functions that create
 * SpillMultiVar objects. These functions use functions that implement variables
 * (return a double) and cuts (return a boolean) to create SpillMultiVar objects
 * that calculate the variables on interactions passing the cuts. The functions
 * are implemented as lambda functions that take a StandardRecord proxy and
 * return a vector of the variable values. Templates are used to handle the
 * different use cases for the variables and cuts.
 * @author mueller@fnal.gov
 */
#ifndef SPINEVAR_H
#define SPINEVAR_H
#include <vector>
#include <iostream>

#include "sbnana/CAFAna/Core/MultiVar.h"

/**
 * @brief Macro to wrap a boolean function in a lambda function.
 * @details This macro is used to wrap a boolean function in a lambda function
 * that it handles the casting of the return value to a double. This is used
 * in the SpineVar functions below to handle the case where the value of a
 * boolean cut is used as a variable.
 */
#define WRAP_BOOL(x) [](auto y) -> double { return x(y); }

#define MCTRUTH caf::Proxy<caf::SRTrueInteraction>
#define TTYPE caf::SRInteractionTruthDLPProxy
#define RTYPE caf::SRInteractionDLPProxy
#define TTYPEP caf::SRParticleTruthDLPProxy
#define RTYPEP caf::SRParticleDLPProxy

/**
 * @brief Function to calculate a variable on each interaction
 * passing some criteria.
 * @details This function is used to calculate a variable on each interaction
 * passing some criteria. The function is implemented as a lambda function that
 * takes a StandardRecord proxy and returns a vector of the variable values.
 * @tparam VARTYPE the type of variable to calculate.
 * @tparam CUTTYPE the type of interaction to loop over.
 * @param fvar the function implementing the variable.
 * @param fcut the function implementing the cut.
 * @param fcat the function implementing the category cut.
 * @return a SpillMultiVar object that calculates the variable.
 */
template<class VARTYPE, class CUTTYPE>
ana::SpillMultiVar SpineVar(double (*fvar)(const VARTYPE &), bool (*fcut)(const CUTTYPE &), bool (*fcat)(const caf::SRInteractionTruthDLPProxy &))
{
    /**
     * @details This is the case that handles the broadcasting of a variable
     * over the SRTrueInteractions matching reco interactions passing a cut.
     * The function iterates over the reco interactions, checks that the
     * interaction passes the cut, that it is matched to a true interaction
     * of a certain truth category, then calls the function implementing the
     * variable on the SRTrueInteraction containing the neutrino information
     * for the true interaction. The result is stored in the vector of doubles.
     */
    if constexpr(std::is_same_v<CUTTYPE, RTYPE> && std::is_same_v<VARTYPE, MCTRUTH>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            bool is_mc(sr->ndlp_true != 0);
            for(auto const& i : sr->dlp)
            {
                if(fcut(i) && i.match_ids.size() > 0 && fcat(sr->dlp_true[i.match_ids[0]]))
                {
                    int64_t nu_id = sr->dlp_true[i.match_ids[0]].nu_id;
                    var.push_back(nu_id >= 0 ? fvar(sr->mc.nu[nu_id]) : PLACEHOLDERVALUE);
                }
                else if(fcut(i) && !is_mc)
                    var.push_back(PLACEHOLDERVALUE);
            }
            return var;
        });
    }

    /**
     * @details This is the case that handles the broadcasting of a variable
     * over the SRTrueInteractions associated with true interactions passing
     * a cut. The function iterates over the true interactions, checks that
     * the interaction passes the cut and belongs to the specified truth
     * category, and that it is matched to a reco interaction. If these
     * conditions are met, the function implementing the variable is called
     * on the SRTrueInteraction containing the neutrino information for the
     * true interaction. The result is stored in the vector of doubles.
     */
    else if constexpr(std::is_same_v<CUTTYPE, TTYPE> && std::is_same_v<VARTYPE, MCTRUTH>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            for(auto const& i : sr->dlp_true)
            {
                if(fcut(i) && i.match_ids.size() > 0 && fcat(i))
                {
                    int64_t nu_id = i.nu_id;
                    var.push_back(nu_id >= 0 ? fvar(sr->mc.nu[nu_id]) : PLACEHOLDERVALUE);
                }
            }
            return var;
        });
    }

    /**
     * @details This is the case that handles the broadcasting of a variable
     * over the true interactions passing a category cut (e.g. the neutrino)
     * that are matched to by reco interactions passing a cut. The function
     * iterates over the reco interactions, checks that the interaction passes
     * the cut, that it is matched to a true interaction, and that the true
     * interaction passes the category cut (or that the record is data). If these
     * conditions are met, the function retrieves the index of the true interaction
     * from the matched reco interaction and calls the function implementing the
     * variable on the true interaction. The result is stored in the vector of
     * variables.
     */
    else if constexpr(std::is_same_v<CUTTYPE, RTYPE> && std::is_same_v<VARTYPE, TTYPE>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            bool is_mc(sr->ndlp_true != 0);
            for(auto const& i : sr->dlp)
            {
                if(fcut(i) && ((i.match_ids.size() > 0 && fcat(sr->dlp_true[i.match_ids[0]])) || !is_mc))
                    var.push_back(i.match_ids.size() > 0 ? fvar(sr->dlp_true[i.match_ids[0]]) : PLACEHOLDERVALUE);
            }
            return var;
        });
    }
    
    /**
     * @details This is the case that handles the broadcasting of a variable
     * over the reco interactions passing a cut that are matched to by true
     * interactions passing a category cut. The function iterates over the reco
     * interactions, checks that the interaction passes the cut, that it is
     * matched to a true interaction, and that the true interaction passes the
     * category cut (or that the record is data). If these conditions are met,
     * the function retrieves the index of the reco interaction from the matched
     * true interaction and calls the function implementing the variable on the
     * reco interaction. The result is stored in the vector of variables.
     */
    else if constexpr(std::is_same_v<CUTTYPE, RTYPE> && std::is_same_v<VARTYPE, RTYPE>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            bool is_mc(sr->ndlp_true != 0);
            for(auto const& i : sr->dlp)
            {
                if(fcut(i) && ((i.match_ids.size() > 0 && fcat(sr->dlp_true[i.match_ids[0]])) || !is_mc))
                    var.push_back(fvar(i));
            }
            return var;
        });
    }

    /**
     * @brief This is the case that handles the broadcasting of a variable
     * over the reco interactions that are matched to by true interactions
     * passing a cut. The function iterates over the true interactions, checks
     * that the interaction passes the cut and that it is matched to a reco
     * interaction. If these conditions are met, the function retrieves the index
     * of the reco interaction from the matched true interaction and calls the
     * function implementing the variable on the reco interaction. The result is
     * stored in the vector of variables.
     */
    else if constexpr(std::is_same_v<CUTTYPE, TTYPE> && std::is_same_v<VARTYPE, RTYPE>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            for(auto const& i : sr->dlp_true)
            {
                if(fcut(i) && i.match_ids.size() > 0)
                    var.push_back(fvar(sr->dlp[i.match_ids[0]]));
            }
            return var;
        });
    }

    /**
     * @brief This is the case that handles the broadcasting of a variable
     * over the true interactions that are matched to by true interactions
     * passing a cut. The function iterates over the true interactions, checks
     * that the interaction passes the cut and that it is matched to a reco
     * interaction. If these conditions are met, the function implementing the
     * variable is called on the true interaction. The result is stored in the
     * vector of variables.
     */
    else if constexpr(std::is_same_v<CUTTYPE, TTYPE> && std::is_same_v<VARTYPE, TTYPE>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            for(auto const& i : sr->dlp_true)
            {
                if(fcut(i) && i.match_ids.size() > 0)
                    var.push_back(fvar(i));
            }
            return var;
        });
    }
    else return ana::SpillMultiVar([](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>{return {1.0};});
}

/**
 * @brief Function to calculate a variable on a single particle in an
 * interaction.
 * @details This function is used to calculate a variable on a single particle
 * in an interaction. The function is templated on the type of particle, the
 * type of interaction to loop over, and the type of particle that is being 
 * that is being used to identify the particle of interest. Practically, this
 * matters because one may wish to identify a true particle and then apply a
 * variable to the corresponding reco particle. The function is implemented as
 * a lambda function that takes a StandardRecord proxy and returns a vector of
 * the variable values.
 * @tparam VARTYPE the type of particle to apply the variable on.
 * @tparam CUTTYPE the type of interaction to loop over.
 * @tparam U the type of particle that is being used to identify the particle of
 * interest.
 * @param fvar the function implementing the variable.
 * @param fcut the function implementing the cut.
 * @param fcat the function implementing the category cut.
 * @param pident the function to identify the particle of interest (by index).
 * @return a SpillMultiVar object that calculates the variable.
 */
template<class VARTYPE, class CUTTYPE, class U>
ana::SpillMultiVar SpineVar(double (*fvar)(const VARTYPE &), bool (*fcut)(const CUTTYPE &), bool (*fcat)(const caf::SRInteractionTruthDLPProxy &), size_t (*pident)(const U &))
{
    /**
     * @details This is the case that handles the mapping of an identified true
     * particle to its corresponding reco particle twin. The function iterates
     * over the reco interactions, checks that the interaction passes the cut,
     * that it is matched to a true interaction, and that the true interaction
     * passes the category cut (or that the record is data). If these conditions
     * are met, the function retrieves the index of the identified true particle
     * from the matched true interaction and calls the function implementing the
     * variable on the true particle. The result is stored in the vector of
     * variables.
     */
    if constexpr(std::is_same_v<CUTTYPE, RTYPE> && std::is_same_v<VARTYPE, TTYPEP>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat, pident](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            bool is_mc(sr->ndlp_true != 0);
            for(auto const& i : sr->dlp)
            {
                if(fcut(i) && i.match_ids.size() > 0 && fcat(sr->dlp_true[i.match_ids[0]]))
                {
                    size_t index(pident(sr->dlp_true[i.match_ids[0]]));
                    var.push_back(i.match_ids.size() > 0 ? fvar(sr->dlp_true[i.match_ids[0]].particles[index]) : PLACEHOLDERVALUE);
                }
                else if(fcut(i) && !is_mc)
                    var.push_back(PLACEHOLDERVALUE);
            }
            return var;
        });
    }
    
    /**
     * @details This is the case that handles the mapping of a reco particle onto
     * some variable. The function iterates over the reco interactions, checks
     * that the interaction passes the cut, that it is matched to a true
     * interaction, and that the true interaction passes the category cut (or 
     * that the record is data). If these conditions are met, the function
     * retrieves the index of the identified particle and calls the function
     * implementing the variable on the identified reco particle. The result is
     * stored in the vector of variables.
     */
    else if constexpr(std::is_same_v<CUTTYPE, RTYPE> && std::is_same_v<VARTYPE, RTYPEP>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat, pident](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            bool is_mc(sr->ndlp_true != 0);
            for(auto const& i : sr->dlp)
            {
                if(fcut(i) && ((i.match_ids.size() > 0 && fcat(sr->dlp_true[i.match_ids[0]])) || !is_mc))
                {
                    size_t index(pident(i));
                    var.push_back(fvar(i.particles[index]));
                }
            }
            return var;
        });
    }
    
    /**
     * @details This is the case that handles the mapping of an identified true
     * particle onto some variable. The function iterates over the true
     * interactions, checks that the interaction passes the cut and that it is
     * matched to a reco interaction. If these conditions are met, the function
     * retrieves the index of the identified particle and calls the function
     * implementing the variable on the identified true particle. The result is
     * stored in the vector of variables.
     */
    else if constexpr(std::is_same_v<CUTTYPE, TTYPE> && std::is_same_v<VARTYPE, TTYPEP>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat, pident](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            for(auto const& i : sr->dlp_true)
            {
                if(fcut(i) && i.match_ids.size() > 0)
                {
                    size_t index(pident(i));
                    var.push_back(fvar(i.particles[index]));
                }
            }
            return var;
        });
    }
    /**
     * @details This is the case that handles the mapping of an identified true
     * particle to its corresponding reco particle twin. It is necessary to
     * first build a map of reco particles by their ID to allow for quick
     * access to the matched reco particle. Then, the function iterates over
     * the true interactions, checks that the interaction passes the cut and
     * that it is matched to a reco interaction. If these conditions are met,
     * the function retrieves the index of the identified particle and grabs
     * the corresponding matched reco particle from the map, if it exists. The
     * function implementing the variable is then called on the reco particle
     * and the result is stored in the vector of variables.
     */
    else if constexpr(std::is_same_v<CUTTYPE, TTYPE> && std::is_same_v<VARTYPE, RTYPEP>)
    {
        return ana::SpillMultiVar([fvar, fcut, fcat, pident](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
        {
            std::vector<double> var;
            std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleDLP> *> reco_particles;
            for(auto const& i : sr->dlp)
            {
                for(auto const& j : i.particles)
                {
                    reco_particles.insert(std::make_pair(j.id, &j));
                }
            }
            for(auto const& i : sr->dlp_true)
            {
                if(fcut(i) && i.match_ids.size() > 0)
                {
                    size_t index(pident(i));
                    if(i.particles[index].match_ids.size() > 0)
                        var.push_back(fvar(*reco_particles[i.particles[index].match_ids[0]]));
                    else var.push_back(PLACEHOLDERVALUE);
                }
            }
            return var;
        });
    }
    else return ana::SpillMultiVar([](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>{return {1.0};});
}

template<class VARTYPE, class PCUTTYPE, class ICUTTYPE>
ana::SpillMultiVar SpineVar(double (*fvar)(const VARTYPE &), bool (*fpcut)(const PCUTTYPE &), bool (*ficut)(const ICUTTYPE &))
{
    /**
     * What do we want to do for p-vars?
     * Select particles (pcut) and the calculate the variable (pvar) for each
     * particle in the interaction.
     * ICUTTYPE tells you which interactions to loop over.
     * PCUTTYPE tells you which particles to select (should match true/reco value with ICUTTYPE)
     * VARTYPE tells you if you apply the variable to the particle itself or the matched one.
     */

    /**
     * @details This is the case that handles the broadcasting of a particle
     * variable over some selection of reco particles 
     */
    if constexpr(std::is_same_v<ICUTTYPE, RTYPE> && std::is_same_v<PCUTTYPE, RTYPEP>)
    {
        if constexpr(std::is_same_v<VARTYPE, RTYPEP>)
        {
            return ana::SpillMultiVar([fvar, fpcut, ficut](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
            {
                std::vector<double> var;
                bool is_mc(sr->ndlp_true != 0);
                for(auto const& i : sr->dlp)
                {
                    if(ficut(i))
                    {
                        for(auto const & j : i.particles)
                        {
                            if(fpcut(j) && (j.match_ids.size() > 0 || !is_mc))
                                var.push_back(fvar(j));
                        }
                    }
                }
                return var;
            });
        }
        else if constexpr(std::is_same_v<VARTYPE, TTYPEP>)
        {
            return ana::SpillMultiVar([fvar, fpcut, ficut](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
            {
                std::vector<double> var;
                std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleTruthDLP> *> true_particles;
                for(auto const& i : sr->dlp_true)
                {
                    for(auto const& j : i.particles)
                    {
                        true_particles.insert(std::make_pair(j.id, &j));
                    }
                }
                for(auto const& i : sr->dlp)
                {
                    if(ficut(i))
                    {
                        for(auto const & j : i.particles)
                        {
                            if(fpcut(j) && j.match_ids.size() > 0)
                                var.push_back(fvar(*true_particles[j.match_ids[0]]));
                        }
                    }
                }
                return var;
            });
        }
    }
    else if constexpr(std::is_same_v<ICUTTYPE, TTYPE> && std::is_same_v<PCUTTYPE, TTYPEP>)
    {
        if constexpr(std::is_same_v<VARTYPE, RTYPEP>)
        {
            return ana::SpillMultiVar([fvar, fpcut, ficut](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
            {
                std::vector<double> var;
                std::map<caf::Proxy<int64_t>, const caf::Proxy<caf::SRParticleDLP> *> reco_particles;
                bool is_mc(sr->ndlp_true != 0);
                for(auto const& i : sr->dlp)
                {
                    for(auto const& j : i.particles)
                    {
                        reco_particles.insert(std::make_pair(j.id, &j));
                    }
                }
                for(auto const& i : sr->dlp_true)
                {
                    if(ficut(i))
                    {
                        for(auto const & j : i.particles)
                        {
                            if(fpcut(j) && (j.match_ids.size() > 0 || !is_mc))
                                var.push_back(fvar(*reco_particles[j.match_ids[0]]));
                            else if(fpcut(j) && is_mc)
                                var.push_back(PLACEHOLDERVALUE);
                        }
                    }
                }
                return var;
            });
        }
        else if constexpr(std::is_same_v<VARTYPE, TTYPEP>)
        {
            return ana::SpillMultiVar([fvar, fpcut, ficut](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>
            {
                std::vector<double> var;
                for(auto const& i : sr->dlp_true)
                {
                    if(ficut(i))
                    {
                        for(auto const & j : i.particles)
                        {
                            if(fpcut(j) && j.match_ids.size() > 0)
                                var.push_back(fvar(j));
                            else if(fpcut(j))
                                var.push_back(fvar(j));
                        }
                    }
                }
                return var;
            });
        }
    }
    else return ana::SpillMultiVar([](const caf::Proxy<caf::StandardRecord> * sr) -> std::vector<double>{return {1.0};});
}
#endif // SPINEVAR_H