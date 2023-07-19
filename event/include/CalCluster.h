/**
 * @file CalCluster.cxx
 * @brief Class used to encapsulate calorimeter cluster information.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#ifndef _CAL_CLUSTER_H__
#define _CAL_CLUSTER_H__

//----------------//
//   C++ StdLib   //
//----------------//
#include <vector>

//----------//
//   ROOT   //
//----------//
#include <TClonesArray.h>
#include <TObject.h>
#include <TRefArray.h>
#include <TRef.h>

class CalCluster : public TObject { 

    public:

        /** Constructor */
        CalCluster();

        /** Destructor */
        ~CalCluster();

        /** Reset the Cluster object */ 
        void Clear(Option_t *option="");

        /**
         * Add a reference to a calorimeter hit composing this cluster.
         *
         * @param hit : Cal hit composing with this cluster
         */
        void addHit(std::size_t i); 

        /** 
         * @return An array of references to the calorimeter hits composing 
         * this cluster. 
         */
        const std::vector<std::size_t>& getHits() const { return hits_; }

        /**
         * @return number of references to the calorimeter hits composing
         * this cluster.
         */
        std::size_t getNHits() const { return hits_.size(); }

        /**
         * Set the position of the calorimeter cluster.
         *
         * @param position : The position of the calorimeter cluster
         */
        void setPosition(const float* position);

        /** @return The position of the calorimeter cluster. */
        std::vector<double> getPosition() const { return { x_, y_, z_ }; };  

        /**
         * Set the energy of the calorimeter cluster.
         *
         * @param energy : The energy of the calorimeter cluster.
         */
        void setEnergy(const double energy) { energy_ = energy; };

        /** @return The energy of the calorimeter cluster. */
        double getEnergy() const { return energy_; };

        /** 
         * Set the time of the calorimeter clusters. 
         *
         * @param time The cluster time
         */
        void setTime(const double time) { time_ = time; }; 

        /** @return The time of the cluster. */
        double getTime() const { return time_; };

        /** 
         * Set the cluster seed i.e. the hit with the highest energy.
         *
         * @param seed The cluster seed. 
         */
        void setSeed(std::size_t i) { seed_hit_ = i; }; 

        /** @return The seed hit of the cluster. */
        std::size_t getSeed() const { return seed_hit_; };

        ClassDef(CalCluster, 1);	

    private:

        /** An array of references to the hits associated withi this cluster. */        
        std::vector<std::size_t> hits_{};

        /** A reference to the seed hit of this cluster. */ 
        std::size_t seed_hit_; 

        /** The x position of the cluster in (mm). */
        double x_{-9999}; 

        /** The y position of the cluster in (mm). */
        double y_{-9999}; 

        /** The z position of the cluster in (mm). */
        double z_{-9999};

        /** The energy of the cluster in GeV. */
        double energy_{-9999}; 

        /** The cluster time. */ 
        double time_{-9999};

}; // CalCluster

#endif // _CALORIMETER_CLUSTER_H_
