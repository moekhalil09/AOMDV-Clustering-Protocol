/*
 * Copyright (c) 2008, Marcello Caleffi, <marcello.caleffi@unina.it>,
 * http://wpage.unina.it/marcello.caleffi
 *
 * The AOMDV code has been developed at DIET, Department of Electronic
 * and Telecommunication Engineering, University of Naples "Federico II"
 *
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License,
 * version 2, as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
 *
 * The copyright of this module includes the following
 * linking-with-specific-other-licenses addition:
 *
 * In addition, as a special exception, the copyright holders of
 * this module give you permission to combine (via static or
															  * dynamic linking) this module with free software programs or
 * libraries that are released under the GNU LGPL and with code
 * included in the standard release of ns-2 under the Apache 2.0
 * license or under otherwise-compatible licenses with advertising
 * requirements (or modified versions of such code, with unchanged
					  * license).  You may copy and distribute such a system following the
 * terms of the GNU GPL for this module and the licenses of the
 * other code concerned, provided that you include the source code of
 * that other code when and as the GNU GPL requires distribution of
 * source code.
 *
 * Note that people who make modified versions of this module
 * are not obligated to grant this special exception for their
 * modified versions; it is their choice whether to do so.  The GNU
 * General Public License gives permission to release a modified
 * version without this exception; this exception also makes it
 * possible to release a modified version which carries forward this
 * exception.
 *
 */



/*
 Copyright (c) 1997, 1998 Carnegie Mellon University.  All Rights
 Reserved. 
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 3. The name of the author may not be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 The AODV code developed by the CMU/MONARCH group was optimized and tuned by Samir Das and Mahesh Marina, University of Cincinnati. The work was partially done in Sun Microsystems.
 
 */



//#include <ip.h>

#include <aomdv/aomdv.h>
#include <aomdv/aomdv_packet.h>
#include <random.h>
#include <cmu-trace.h>
#include <vector>
#include <cstdio>
#include <numeric>
//#include <energy-model.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <limits>
#include <cstdlib>
#include <algorithm>
#define max(a,b)        ( (a) > (b) ? (a) : (b) )
#define CURRENT_TIME    Scheduler::instance().clock()
using namespace std;
// AOMDV
#define DROP_RTR_RTEXPIRE               "REXP"
#define DROP_RTR_HELLO                  "HELO"

//#define DEBUG
//#define ERROR

#ifdef DEBUG
//static int extra_route_reply = 0;
//static int limit_route_request = 0;
static int route_request = 0;
#endif


/*
 TCL Hookss
 */

const int numNodes = 150;
std::vector<float> energy(numNodes, 5000); 
std::vector<float> current_time(numNodes, 180); 
std::vector<bool> changed(numNodes,false);
float sum;
float mean;
float sumcr;
float meancr;
int off_nodes = 0;
bool itr = true;
int cpt = 0;
std::vector<int> CH;
std::vector<int> clusters(numNodes, -1);
std::vector<int> final_clusters(numNodes+1, -1);
std::vector<int> final_CH;
std::vector<int> dead_ch;
float energy_consome = 0;
float dv = 180.00;
bool once = true;
std::vector<std::pair<double, double> > energy_consome_s;

/* ENERGY MODEL */
// Constants
const float E_elec = 50e-4;  // 50 nJ/bit
const float E_amp = 1e-12;  // 100 pJ/bit/m^2
const int RREQ_SIZE_BITS = 240;  // 30 bytes = 240 bits
const int RREP_SIZE_BITS = 240;  // 30 bytes = 240 bits
const int SEND_REQ_SIZE_BITS = 400;  // 50 bytes = 400 bits
const int SEND_REP_SIZE_BITS = 400;  // 50 bytes = 400 bits
float rayon = 20.0;  // Distance in meters

// Calculate transmission energy
float calculateTxEnergy(int messageSizeBits, float distance) {
    return (E_elec * messageSizeBits) + (E_amp * messageSizeBits * distance * distance);
}

// Calculate reception energy
float calculateRxEnergy(int messageSizeBits) {
    return E_elec * messageSizeBits;
}



/* CLUSTERING FUNCTIONS HANDELING BY MOHAMMED KHALIL*/
/*DISTANCE CALCULATION */
double calculateDistance(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        sum += pow(a[i] - b[i], 2);
    }
    return sqrt(sum);
}
void assignToNearestClusterHead(const std::vector<std::vector<double> >& data, const std::vector<int>& CH, std::vector<int>& clusters) 
 {
    for (size_t i = 0; i < data.size(); ++i) {
        double minDist =  99999999999;
        int nearestClusterHeadIndex = -1;

        // Find the nearest cluster head for the current data point
        for (size_t j = 0; j < CH.size(); ++j) {
            double dist = calculateDistance(data[i], data[CH[j]]);
            if (dist < minDist) {
                minDist = dist;
                nearestClusterHeadIndex = CH[j];
            }
        }

        // Assign the data point to the nearest cluster head
        clusters[i] = nearestClusterHeadIndex;
    }
}
void clusterNodes(const std::vector<std::vector<double> >& distances, double K, std::vector<int>& CH) {
    int n = distances.size(); // Size of the distance matrix (assuming square matrix)

    // Always add nodes 0 and 1 initially to the cluster heads (CH)
    CH.push_back(0);
    CH.push_back(1);
    final_CH.push_back(1);
    final_CH.push_back(2);

    // Loop through the nodes (assuming node IDs are indices)
    for (int j = 2; j < n; ++j) {
        bool isClusterHead = true; // Assume j is a new cluster head

        // Check distances from node j to existing cluster heads in CH
        for (size_t k = 0; k < CH.size(); ++k) {
            int ch = CH[k];
            if (ch < j && distances[j][ch] <= K) {
                isClusterHead = false; // j is not a new cluster head
                break;
            }
        }

        // If j is a new cluster head, add it to CH
        if (isClusterHead) {

            CH.push_back(j);
            final_CH.push_back(j+1);

        }
    }

    std::cout << CH.size() << " clusterheads generated." << std::endl;
}
bool is_member_of(int nodeIndex, int clusterHead) {
    return final_clusters[nodeIndex] == clusterHead;
}

bool is_cluster_head(int nodeIndex, const std::vector<int>& CH) {
    for (std::vector<int>::const_iterator it = CH.begin(); it != CH.end(); ++it) {
        if (nodeIndex == *it) {
            return true; //
            std::cout<<"Node "<< nodeIndex<<" is a cluster head";
        }
    }
    return false; std::cout<<"Node "<< nodeIndex<<" is not a cluster head";
}
/*HERE IS THE PART OF THE RECONFIGURATION WHEN A CLUSTERHEAD DIES*/
void reconfiguration(int cluster_head) {
    std::cout << "Reconfiguration of the node " << cluster_head<<std::endl;
    std::vector<int> temp_ch;
    float temp_max = -999999.0f; // Initialize with a float value
    int temp_index = -1; // Initialize with an invalid index

    // Create a temporary list of cluster members belonging to the cluster head
    //std::cout<<"{";
    for (std::vector<int>::size_type i = 0; i < final_clusters.size(); ++i) {
        if (final_clusters[i] == cluster_head) {
            temp_ch.push_back(i);
            //std::cout<<i<<",";
        }
    }
    //std::cout<<"} ";

    // Loop through the temporary list to find the node with maximum energy
    for (std::vector<int>::size_type node = 0; node < temp_ch.size(); ++node) {
    if (energy[temp_ch[node]] > temp_max && static_cast<int>(node) != cluster_head) {
        temp_max = energy[temp_ch[node]];
        temp_index = temp_ch[node]; 
       // std::cout << "node " << temp_ch[node] << " has " << energy[temp_ch[node]] << " energy ";
    }
     
}
    std::cout<<std::endl;
    if(temp_index != -1)
{
    std::cout << temp_index << " will replace the CH: " <<cluster_head<< std::endl;

    // Update the CH list and the clusters
    for (std::vector<int>::size_type i = 0; i < final_CH.size(); ++i) {
        if (final_CH[i] == cluster_head) {
            final_CH[i] = temp_index;
            break; // Assuming there's only one occurrence of cluster_head in CH
        }
    }
    
    // Update the clusters
    for (std::vector<int>::size_type i = 0; i < final_clusters.size(); ++i) {
        if (final_clusters[i] == cluster_head) {
            final_clusters[i] = temp_index;
        }
    }
    //update the dead cluster list
    if(temp_index != cluster_head){
    	dead_ch.erase(std::remove(dead_ch.begin(), dead_ch.end(), cluster_head), dead_ch.end());
    }
}

}


int hdr_aomdv::offset_;
static class AOMDVHeaderClass : public PacketHeaderClass {
public:
	AOMDVHeaderClass() : PacketHeaderClass("PacketHeader/AOMDV",
														sizeof(hdr_all_aomdv)) {
		bind_offset(&hdr_aomdv::offset_);
   } 
} class_rtProtoAOMDV_hdr;

static class AOMDVclass : public TclClass {
public:
	AOMDVclass() : TclClass("Agent/AOMDV") {}
	TclObject* create(int argc, const char*const* argv) {
		assert(argc == 5);
		//return (new AODV((nsaddr_t) atoi(argv[4])));
		return (new AOMDV((nsaddr_t) Address::instance().str2addr(argv[4])));
	}
} class_rtProtoAOMDV;


int
AOMDV::command(int argc, const char*const* argv) {
	if(argc == 2) {
		Tcl& tcl = Tcl::instance();
		
		if(strncasecmp(argv[1], "id", 2) == 0) {
			tcl.resultf("%d", index);
			return TCL_OK;
		}
		// AOMDV code - should it be removed?
		if (strncasecmp(argv[1], "dump-table", 10) == 0) {
			printf("Node %d: Route table:\n", index);
			rtable.rt_dumptable();
			return TCL_OK;
		}    
		if(strncasecmp(argv[1], "start", 2) == 0) {
			btimer.handle((Event*) 0);
			
#ifndef AOMDV_LINK_LAYER_DETECTION
			htimer.handle((Event*) 0);
			ntimer.handle((Event*) 0);
#endif // LINK LAYER DETECTION
			
			rtimer.handle((Event*) 0);
			return TCL_OK;
		}
	}
	else if(argc == 3) {
		if(strcmp(argv[1], "index") == 0) {
			index = atoi(argv[2]);
			return TCL_OK;
		}
		
		else if(strcmp(argv[1], "log-target") == 0 || strcmp(argv[1], "tracetarget") == 0) {
			logtarget = (Trace*) TclObject::lookup(argv[2]);
			if(logtarget == 0)
				return TCL_ERROR;
			return TCL_OK;
		}
		else if(strcmp(argv[1], "drop-target") == 0) {
			int stat = rqueue.command(argc,argv);
			if (stat != TCL_OK) return stat;
			return Agent::command(argc, argv);
		}
		else if(strcmp(argv[1], "if-queue") == 0) {
			AOMDVifqueue = (PriQueue*) TclObject::lookup(argv[2]);
			
			if(AOMDVifqueue == 0)
				return TCL_ERROR;
			return TCL_OK;
		}
		// AODV ns-2.31 code
		else if (strcmp(argv[1], "port-dmux") == 0) {
			dmux_ = (PortClassifier *)TclObject::lookup(argv[2]);
			if (dmux_ == 0) {
				fprintf (stderr, "%s: %s lookup of %s failed\n", __FILE__,
							argv[1], argv[2]);
				return TCL_ERROR;
			}
			return TCL_OK;
		}
	}
	return Agent::command(argc, argv);
}

/* 
Constructor
 */

AOMDV::AOMDV(nsaddr_t id) : Agent(PT_AOMDV),
btimer(this), htimer(this), ntimer(this), 
rtimer(this), lrtimer(this), rqueue() {
	
	// AOMDV code
	aomdv_max_paths_ = 3;
	bind("aomdv_max_paths_", &aomdv_max_paths_);
	
	aomdv_prim_alt_path_len_diff_ = 1;
	bind("aomdv_prim_alt_path_len_diff_", &aomdv_prim_alt_path_len_diff_);
	
	index = id;
	seqno = 2;
	bid = 1;
	
	LIST_INIT(&nbhead);
	LIST_INIT(&bihead);
	
	logtarget = 0;
	AOMDVifqueue = 0;
}

/*
 Timers
 */

void
AOMDVBroadcastTimer::handle(Event*) {
	agent->id_purge();
	Scheduler::instance().schedule(this, &intr, BCAST_ID_SAVE);
}

void
AOMDVHelloTimer::handle(Event*) {
	// AOMDV code - could it be removed?
	// agent->sendHello();
	/* Do not send a HELLO message unless we have a valid route entry. */
	if (agent->rtable.rt_has_active_route())
		agent->sendHello();
	// AODV ns-2.31 code
	// double interval = HELLO_INTERVAL + 0.01*Random::uniform();
	double interval = MinHelloInterval + ((MaxHelloInterval - MinHelloInterval) * Random::uniform());
	assert(interval >= 0);
	Scheduler::instance().schedule(this, &intr, interval);
}

void
AOMDVNeighborTimer::handle(Event*) {
	agent->nb_purge();
	Scheduler::instance().schedule(this, &intr, HELLO_INTERVAL);
}

void
AOMDVRouteCacheTimer::handle(Event*) {
	agent->rt_purge();	
#define FREQUENCY 0.5 // sec
	Scheduler::instance().schedule(this, &intr, FREQUENCY);
}

void
AOMDVLocalRepairTimer::handle(Event* p)  {  // SRD: 5/4/99
	aomdv_rt_entry *rt;
	struct hdr_ip *ih = HDR_IP( (Packet *)p);
	
	/* you get here after the timeout in a local repair attempt */
	/* fprintf(stderr, "%s\n", __FUNCTION__); */
	
	
	rt = agent->rtable.rt_lookup(ih->daddr());
	
	if (rt && rt->rt_flags != RTF_UP) {
		// route is yet to be repaired
		// I will be conservative and bring down the route
		// and send route errors upstream.
		/* The following assert fails, not sure why */
		/* assert (rt->rt_flags == RTF_IN_REPAIR); */
		
		//rt->rt_seqno++;
		agent->rt_down(rt);
		// send RERR
#ifdef DEBUG
//		fprintf(stderr,"Node %d: Dst - %d, failed local repair\n",index, rt->rt_dst);
#endif      
	}
	Packet::free((Packet *)p);
}


/*
 Broadcast ID Management  Functions
 */


// AODV ns-2.31 code
void
AOMDV::id_insert(nsaddr_t id, u_int32_t bid) {
	AOMDVBroadcastID *b = new AOMDVBroadcastID(id, bid);
	
	assert(b);
	b->expire = CURRENT_TIME + BCAST_ID_SAVE;
	LIST_INSERT_HEAD(&bihead, b, link);
}

// AODV ns-2.31 code
/* SRD */
bool
AOMDV::id_lookup(nsaddr_t id, u_int32_t bid) {
	AOMDVBroadcastID *b = bihead.lh_first;
	
	// Search the list for a match of source and bid
	for( ; b; b = b->link.le_next) {
		if ((b->src == id) && (b->id == bid))
			return true;     
	}
	return false;
}

// AOMDV ns-2.31 code
AOMDVBroadcastID*
AOMDV::id_get(nsaddr_t id, u_int32_t bid) {
	AOMDVBroadcastID *b = bihead.lh_first;
	
	// Search the list for a match of source and bid
	for( ; b; b = b->link.le_next) {
		if ((b->src == id) && (b->id == bid))
			return b;     
	}
	return NULL;
}

void
AOMDV::id_purge() {
	AOMDVBroadcastID *b = bihead.lh_first;
	AOMDVBroadcastID *bn;
	double now = CURRENT_TIME;
	
	for(; b; b = bn) {
		bn = b->link.le_next;
		if(b->expire <= now) {
			LIST_REMOVE(b,link);
			delete b;
		}
	}
}

/*
 Helper Functions
 */

double
AOMDV::PerHopTime(aomdv_rt_entry *rt) {
	int num_non_zero = 0, i;
	double total_latency = 0.0;
	
	if (!rt)
		return ((double) NODE_TRAVERSAL_TIME );
	
	for (i=0; i < MAX_HISTORY; i++) {
		if (rt->rt_disc_latency[i] > 0.0) {
			num_non_zero++;
			total_latency += rt->rt_disc_latency[i];
		}
	}
	if (num_non_zero > 0)
		return(total_latency / (double) num_non_zero);
	else
		return((double) NODE_TRAVERSAL_TIME);
	
}

/*
 Link Failure Management Functions
 */

static void
aomdv_rt_failed_callback(Packet *p, void *arg) {
	((AOMDV*) arg)->rt_ll_failed(p);
}

/*
 * This routine is invoked when the link-layer reports a route failed.
 */
void
AOMDV::rt_ll_failed(Packet *p) {
// AOMDV ns-2.31 code
#ifndef AOMDV_LINK_LAYER_DETECTION
	drop(p, DROP_RTR_MAC_CALLBACK);
#else 
	
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	aomdv_rt_entry *rt;
	nsaddr_t broken_nbr = ch->next_hop_;
	
	/*
	 * Non-data packets and Broadcast Packets can be dropped.
	 */
	if(! DATA_PACKET(ch->ptype()) ||
		(u_int32_t) ih->daddr() == IP_BROADCAST) {
		drop(p, DROP_RTR_MAC_CALLBACK);
		return;
	}
	log_link_broke(p);
	if((rt = rtable.rt_lookup(ih->daddr())) == 0) {
		drop(p, DROP_RTR_MAC_CALLBACK);
		return;
	}
	log_link_del(ch->next_hop_);
	
#ifdef AOMDV_LOCAL_REPAIR
	/* if the broken link is closer to the dest than source, 
		attempt a local repair. Otherwise, bring down the route. */
	
	
	if (ch->num_forwards() > rt->rt_hops) {
		local_rt_repair(rt, p); // local repair
										// retrieve all the packets in the ifq using this link,
										// queue the packets for which local repair is done, 
		return;
	} 
#endif // LOCAL REPAIR  
	
	{
		// AOMDV code - could it be removed?
		handle_link_failure(broken_nbr);
		// AOMDV code
#ifdef AOMDV_PACKET_SALVAGING
		
		if ( !DATA_PACKET(ch->ptype()) ) 
			drop(p, DROP_RTR_MAC_CALLBACK);
		else {
			// salvage the packet using an alternate path if available.
			aomdv_rt_entry *rt = rtable.rt_lookup(ih->daddr());
			if ( rt && (rt->rt_flags == RTF_UP) && (ch->aomdv_salvage_count_ < AOMDV_MAX_SALVAGE_COUNT) ) {
				ch->aomdv_salvage_count_ += 1;
				forward(rt, p, NO_AOMDV_DELAY);
			}
			else drop(p, DROP_RTR_MAC_CALLBACK);
		}
		while((p = AOMDVifqueue->filter(broken_nbr))) {
			struct hdr_cmn *ch = HDR_CMN(p);
			struct hdr_ip *ih = HDR_IP(p);
			if ( !DATA_PACKET(ch->ptype()) ) 
				drop(p, DROP_RTR_MAC_CALLBACK);
			else {
				// salvage the packet using an alternate path if available.
				aomdv_rt_entry *rt = rtable.rt_lookup(ih->daddr());
				if ( rt && (rt->rt_flags == RTF_UP) && (ch->aomdv_salvage_count_ < AOMDV_MAX_SALVAGE_COUNT) ) {
					ch->aomdv_salvage_count_ += 1;
					forward(rt, p, NO_AOMDV_DELAY);
				}
				else drop(p, DROP_RTR_MAC_CALLBACK);
			}
		} 
#else // NO PACKET SALVAGING
		drop(p, DROP_RTR_MAC_CALLBACK);
		// Do the same thing for other packets in the interface queue using the
		// broken link -Mahesh
		while((p = AOMDVifqueue->filter(broken_nbr))) {
			drop(p, DROP_RTR_MAC_CALLBACK);
		} 
		nb_delete(broken_nbr);
		// AOMDV code
#endif // NO PACKET SALVAGING        
	}
	
#endif // LINK LAYER DETECTION
}

// AOMDV code
void
AOMDV::handle_link_failure(nsaddr_t id) {
	bool error=true;
	aomdv_rt_entry *rt, *rtn;
	Packet *rerr = Packet::alloc();
	struct hdr_aomdv_error *re = HDR_AOMDV_ERROR(rerr);
#ifdef DEBUG
	fprintf(stderr, "%s: multipath version\n", __FUNCTION__);
#endif // DEBUG
	re->DestCount = 0;
	for(rt = rtable.head(); rt; rt = rtn) {  // for each rt entry
		AOMDV_Path* path;
		rtn = rt->rt_link.le_next; 
		if ((rt->rt_flags == RTF_UP) && (path=rt->path_lookup(id)) ) {
			assert((rt->rt_seqno%2) == 0);
			
			rt->path_delete(id);
			if (rt->path_empty()) {
				rt->rt_seqno++;
				rt->rt_seqno = max(rt->rt_seqno, rt->rt_highest_seqno_heard);
				// CHANGE
				if (rt->rt_error) {
					re->unreachable_dst[re->DestCount] = rt->rt_dst;
					re->unreachable_dst_seqno[re->DestCount] = rt->rt_seqno;
#ifdef DEBUG
					fprintf(stderr, "%s(%f): %d\t(%d\t%u\t%d)\n", __FUNCTION__, CURRENT_TIME,
							  index, re->unreachable_dst[re->DestCount],
							  re->unreachable_dst_seqno[re->DestCount], id);
#endif // DEBUG
					re->DestCount += 1;
					rt->rt_error = false;
				}
				// CHANGE
				rt_down(rt);
			}
		}
	}   
	
	if ( (re->DestCount > 0) && (error) ) {
#ifdef DEBUG
		fprintf(stdout, "%s(%f): %d\tsending RERR...\n", __FUNCTION__, CURRENT_TIME, index);
#endif // DEBUG
		sendError(rerr, false);
	}
	else {
		Packet::free(rerr);
	}
}

void
AOMDV::local_rt_repair(aomdv_rt_entry *rt, Packet *p) {
#ifdef DEBUG
	fprintf(stderr,"%s: Dst - %d\n", __FUNCTION__, rt->rt_dst); 
#endif  
	// Buffer the packet 
	rqueue.enque(p);
	
	// mark the route as under repair 
	rt->rt_flags = RTF_IN_REPAIR;
	
	sendRequest(rt->rt_dst);
	
	// set up a timer interrupt
	Scheduler::instance().schedule(&lrtimer, p->copy(), rt->rt_req_timeout);
}

void
AOMDV::rt_down(aomdv_rt_entry *rt) {
	/*
	 *  Make sure that you don't "down" a route more than once.
	 */
	
	// AOMDV code
#ifdef DEBUG
	fprintf(stderr, "%s: multipath version\n", __FUNCTION__);
#endif // DEBUG
	
	if(rt->rt_flags == RTF_DOWN) {
		return;
	}
	
	// AODV ns-2.31 code
	// assert (rt->rt_seqno%2); // is the seqno odd?
	// AOMDV code
	rt->rt_flags = RTF_DOWN;
	rt->rt_advertised_hops = INFINITY;
	rt->path_delete();
	rt->rt_expire = 0;
	
} /* rt_down function */

/*
 Route Handling Functions
 */

void
AOMDV::rt_resolve(Packet *p) {
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	aomdv_rt_entry *rt;
	
	/*
	 *  Set the transmit failure callback.  That
	 *  won't change.
	 */
	ch->xmit_failure_ = aomdv_rt_failed_callback;
	ch->xmit_failure_data_ = (void*) this;
   rt = rtable.rt_lookup(ih->daddr());
	if(rt == 0) {
		rt = rtable.rt_add(ih->daddr());
	}
	
	/*
	 * If the route is up, forward the packet 
	 */
   
	if(rt->rt_flags == RTF_UP) {
		assert(rt->rt_hops != INFINITY2);
		forward(rt, p, NO_AOMDV_DELAY);
	}
	/*
	 *  if I am the source of the packet, then do a Route Request.
	 */
   else if(ih->saddr() == index) {
		rqueue.enque(p);
		sendRequest(rt->rt_dst);
	}
	/*
	 *   A local repair is in progress. Buffer the packet. 
	 */
	else if (rt->rt_flags == RTF_IN_REPAIR) {
		rqueue.enque(p);
	}
	
	/*
	 * I am trying to forward a packet for someone else to which
	 * I don't have a route.
	 */
	else {
		Packet *rerr = Packet::alloc();
		struct hdr_aomdv_error *re = HDR_AOMDV_ERROR(rerr);
		/* 
			* For now, drop the packet and send error upstream.
		 * Now the route errors are broadcast to upstream
		 * neighbors - Mahesh 09/11/99
		 */  
		
		assert (rt->rt_flags == RTF_DOWN);
		re->DestCount = 0;
		re->unreachable_dst[re->DestCount] = rt->rt_dst;
		re->unreachable_dst_seqno[re->DestCount] = rt->rt_seqno;
		re->DestCount += 1;
#ifdef DEBUG
		fprintf(stderr, "%s: sending RERR...\n", __FUNCTION__);
#endif
		sendError(rerr, false);
		
		drop(p, DROP_RTR_NO_ROUTE);
	}
	
}

void
AOMDV::rt_purge() {
	aomdv_rt_entry *rt, *rtn;
	double now = CURRENT_TIME;
	double delay = 0.0;
	Packet *p;
	
	for(rt = rtable.head(); rt; rt = rtn) {  // for each rt entry
		rtn = rt->rt_link.le_next;
		// AOMDV code
		// MODIFIED BY US! Added '&& rt-> ...' in if-statement
		if (rt->rt_flags == RTF_UP && (rt->rt_expire < now)) {
			rt->path_purge();
			if (rt->path_empty()) {
				while((p = rqueue.deque(rt->rt_dst))) {
               drop(p, DROP_RTR_RTEXPIRE);
				}
				rt->rt_seqno++;
				rt->rt_seqno = max(rt->rt_seqno, rt->rt_highest_seqno_heard);
				if (rt->rt_seqno%2 == 0) rt->rt_seqno += 1;
				// assert (rt->rt_seqno%2);
				rt_down(rt);
			}
		}
		else if (rt->rt_flags == RTF_UP) {
			// If the route is not expired,
			// and there are packets in the sendbuffer waiting,
			// forward them. This should not be needed, but this extra 
			// check does no harm.
			assert(rt->rt_hops != INFINITY2);
			while((p = rqueue.deque(rt->rt_dst))) {
				forward (rt, p, delay);
				delay += ARP_DELAY;
			}
		} 
		else if (rqueue.find(rt->rt_dst))
			// If the route is down and 
			// if there is a packet for this destination waiting in
			// the sendbuffer, then send out route request. sendRequest
			// will check whether it is time to really send out request
			// or not.
			// This may not be crucial to do it here, as each generated 
			// packet will do a sendRequest anyway.
			
			sendRequest(rt->rt_dst); 
   }
	
}

/*
 Packet Reception Routines
 */

void
AOMDV::recv(Packet *p, Handler*) {
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	
	assert(initialized());
	//assert(p->incoming == 0);
	// XXXXX NOTE: use of incoming flag has been depracated; In order to track direction of pkt flow, direction_ in hdr_cmn is used instead. see packet.h for details.
	
	if(ch->ptype() == PT_AOMDV) {
		ih->ttl_ -= 1;
		recvAOMDV(p);
		return;
	}
	
	
	/*
	 *  Must be a packet I'm originating...
	 */
	if((ih->saddr() == index) && (ch->num_forwards() == 0)) {
		/*
		 * Add the IP Header
		 */
		ch->size() += IP_HDR_LEN;
		// AOMDV code
		ch->aomdv_salvage_count_ = 0;
		// Added by Parag Dadhania && John Novatnack to handle broadcasting
		if ( (u_int32_t)ih->daddr() != IP_BROADCAST)
			ih->ttl_ = NETWORK_DIAMETER;
	}
	/*
	 *  I received a packet that I sent.  Probably
	 *  a routing loop.
	 */
	else if(ih->saddr() == index) {
		drop(p, DROP_RTR_ROUTE_LOOP);
		return;
	}
	/*
	 *  Packet I'm forwarding...
	 */
	else {
		/*
		 *  Check the TTL.  If it is zero, then discard.
		 */
		if(--ih->ttl_ == 0) {
			drop(p, DROP_RTR_TTL);
			return;
		}
	}
	// Added by Parag Dadhania && John Novatnack to handle broadcasting
	if ( (u_int32_t)ih->daddr() != IP_BROADCAST)
		rt_resolve(p);
	else
		forward((aomdv_rt_entry*) 0, p, NO_AOMDV_DELAY);
}


void
AOMDV::recvAOMDV(Packet *p) {
	struct hdr_aomdv *ah = HDR_AOMDV(p);
	// AODV ns-2.31 code
	// struct hdr_ip *ih = HDR_IP(p);
	assert(HDR_IP (p)->sport() == RT_PORT);
	assert(HDR_IP (p)->dport() == RT_PORT);
	
	/*
	 * Incoming Packets.
	 */
	switch(ah->ah_type) {
		
		case AOMDVTYPE_RREQ:
			recvRequest(p);
			break;
			
		case AOMDVTYPE_RREP:
			recvReply(p);
			break;
			
		case AOMDVTYPE_RERR:
			recvError(p);
			break;
			
		case AOMDVTYPE_HELLO:
			recvHello(p);
			break;
			
		default:
			fprintf(stderr, "Invalid AOMDV type (%x)\n", ah->ah_type);
			exit(1);
	}
	
}


// AOMDV

void
AOMDV::recvRequest(Packet *p) {
	struct hdr_ip *ih = HDR_IP(p);
	struct hdr_aomdv_request *rq = HDR_AOMDV_REQUEST(p);
	aomdv_rt_entry *rt;
	AOMDVBroadcastID* b = NULL;
	bool kill_request_propagation = false;
	AOMDV_Path* reverse_path = NULL;
	//int node_index = index + 1;
	//int src_index = rq->rq_src + 1;
	/*
	 * Drop if:
	 *      - I'm the source
	 *      - I recently heard this request.
	 */
	/*
	std::cout << "index: " << index
              << " is CH: " << std::boolalpha << is_cluster_head(index, final_CH)
              << " source is: " << rq->rq_src
              << " member of: " << final_clusters[rq->rq_src] << std::endl;
*/
	

if ((index == 0 && CURRENT_TIME < 12) || // If itr is true, execute
    ((energy[index] > 250 ) &&
    ((index == 0 && is_cluster_head(rq->rq_src, final_CH)) ||
    (is_cluster_head(index, final_CH) && is_cluster_head(rq->rq_src, final_CH)) ||
    (is_cluster_head(index, final_CH) && is_member_of(rq->rq_src, index)))))
    
    {

	if(rq->rq_src == index) {
#ifdef DEBUG
		fprintf(stderr, "%s: got my own REQUEST\n", __FUNCTION__);
#endif // DEBUG
		Packet::free(p);
		return;
	} 
	if((index != 0) && CURRENT_TIME > 11 && dead_ch.size() < final_CH.size()){
		float rxEnergy = calculateRxEnergy(RREQ_SIZE_BITS);
        energy_consome += rxEnergy;
        energy[index] -= rxEnergy;

        float txEnergy = calculateTxEnergy(RREP_SIZE_BITS, rayon);
        energy_consome += txEnergy;
        energy[index] -= txEnergy;
		
	}
	energy_consome_s.push_back(std::make_pair(energy_consome, CURRENT_TIME));
	
   /* If RREQ has already been received - drop it, else remember "RREQ id" <src IP, bcast ID>. */
   if ( (b = id_get(rq->rq_src, rq->rq_bcast_id)) == NULL)  {
		// Cache the broadcast ID
		id_insert(rq->rq_src, rq->rq_bcast_id);
		b = id_get(rq->rq_src, rq->rq_bcast_id);
   }
   else 
		kill_request_propagation = true;
	
   /* If I am a neighbor to the RREQ source, make myself first hop on path from source to dest. */
   if (rq->rq_hop_count == 0) 
		rq->rq_first_hop = index;
	
	/* 
		* We are either going to forward the REQUEST or generate a
	 * REPLY. Before we do anything, we make sure that the REVERSE
	 * route is in the route table.
	 */
	aomdv_rt_entry *rt0; // rt0 is the reverse route 
   
   rt0 = rtable.rt_lookup(rq->rq_src);
   if(rt0 == 0) { /* if not in the route table */
		// create an entry for the reverse route.
		rt0 = rtable.rt_add(rq->rq_src);
   }

	/*
	 * Create/update reverse path (i.e. path back to RREQ source)
	 * If RREQ contains more recent seq number than route table entry - update route entry to source.
	 */
	if (rt0->rt_seqno < rq->rq_src_seqno) {
		rt0->rt_seqno = rq->rq_src_seqno;
		rt0->rt_advertised_hops = INFINITY;
		rt0->path_delete(); // Delete all previous paths to RREQ source 
		rt0->rt_flags = RTF_UP;
		/* Insert new path for route entry to source of RREQ. 
			(src addr, hop count + 1, lifetime, last hop (first hop for RREQ)) */
		reverse_path = rt0->path_insert(ih->saddr(), rq->rq_hop_count+1, CURRENT_TIME + REV_ROUTE_LIFE, rq->rq_first_hop);
		// CHANGE
		rt0->rt_last_hop_count = rt0->path_get_max_hopcount();
		// CHANGE
	}
	/* If a new path with smaller hop count is received 
	(same seqno, better hop count) - try to insert new path in route table. */
	else if ( (rt0->rt_seqno == rq->rq_src_seqno) &&
				 (rt0->rt_advertised_hops > rq->rq_hop_count)
				 ) {
		AOMDV_Path* erp=NULL;
		
		assert(rt0->rt_flags == RTF_UP);  // Make sure path is up
	
		/*
		 * If path already exists - adjust the lifetime of the path.
		 */
		if ((reverse_path = rt0->disjoint_path_lookup(ih->saddr(), rq->rq_first_hop))) {
			assert(reverse_path->hopcount == (rq->rq_hop_count+1));
			reverse_path->expire = max(reverse_path->expire, (CURRENT_TIME + REV_ROUTE_LIFE)); 
		}
		/*
		 * Got a new alternate disjoint reverse path - so insert it.
		 * I.e. no path exists which has RREQ source as next hop and no 
		 * path with RREQ first hop as last hop exists for this route entry.
		 * Simply stated: no path with the same last hop exists already.
		 */
		else if (rt0->new_disjoint_path(ih->saddr(), rq->rq_first_hop)) {
			/* Only insert new path if not too many paths exists for this destination 
			and new path does not differ too much in length compared to previous paths */
			if ( (rt0->rt_num_paths_ < aomdv_max_paths_) &&
				  (((rq->rq_hop_count + 1) - rt0->path_get_min_hopcount()) <= aomdv_prim_alt_path_len_diff_)
				  ) {
				/* Insert new (disjoint) reverse path */
				reverse_path = rt0->path_insert(ih->saddr(), rq->rq_hop_count+1, CURRENT_TIME + REV_ROUTE_LIFE, rq->rq_first_hop);
				// CHANGE
				rt0->rt_last_hop_count = rt0->path_get_max_hopcount();
				// CHANGE
			}
			/* If new path differs too much in length compared to previous paths - drop packet. */
			if (((rq->rq_hop_count + 1) - rt0->path_get_min_hopcount()) > aomdv_prim_alt_path_len_diff_) {
				Packet::free(p);
				return;
			}
		}
		/* (RREQ was intended for me) AND 
			((Path with RREQ first hop as last hop does not exist) OR 
			 (The path exists and has less hop count than RREQ)) - drop packet. 
			Don't know what this case is for... */
		else if ( (rq->rq_dst == index) && 
					 ( ((erp = rt0->path_lookup_lasthop(rq->rq_first_hop)) == NULL) ||
						((rq->rq_hop_count+1) > erp->hopcount)
						)
					 )  {
			Packet::free(p);
			return;
		}
	}
	/* Older seqno (or same seqno with higher hopcount), i.e. I have a 
	more recent route entry - so drop packet.*/
	else {
		Packet::free(p);
		return;
	}

	/* If route is up */   
	if (rt0->rt_flags == RTF_UP) {
		// Reset the soft state 
		rt0->rt_req_timeout = 0.0; 
		rt0->rt_req_last_ttl = 0;
		rt0->rt_req_cnt = 0;
		
		/* 
			* Find out whether any buffered packet can benefit from the 
		 * reverse route.
		 */
		Packet *buffered_pkt;
		while ((buffered_pkt = rqueue.deque(rt0->rt_dst))) {
			if (rt0 && (rt0->rt_flags == RTF_UP)) {
				forward(rt0, buffered_pkt, NO_AOMDV_DELAY);
			}
		}
	}
	/* Check route entry for RREQ destination */
	rt = rtable.rt_lookup(rq->rq_dst);

	/* I am the intended receiver of the RREQ - so send a RREP */ 
	if (rq->rq_dst == index) {
		
		if (seqno < rq->rq_dst_seqno) {
			//seqno = max(seqno, rq->rq_dst_seqno)+1;
			seqno = rq->rq_dst_seqno + 1; //CHANGE (replaced above line with this one)
		}
		/* Make sure seq number is even (why?) */
		if (seqno%2) 
			seqno++;
		
		
		sendReply(rq->rq_src,              // IP Destination
					 0,                       // Hop Count
					 index,                   // (RREQ) Dest IP Address 
					 seqno,                   // Dest Sequence Num
					 MY_ROUTE_TIMEOUT,        // Lifetime
					 rq->rq_timestamp,        // timestamp
					 ih->saddr(),             // nexthop
					 rq->rq_bcast_id,         // broadcast id to identify this route discovery
					 ih->saddr());         
		
		Packet::free(p);
	}
	/* I have a fresh route entry for RREQ destination - so send RREP */
	else if ( rt &&
				 (rt->rt_flags == RTF_UP) &&
				 (rt->rt_seqno >= rq->rq_dst_seqno) ) {
		
		assert ((rt->rt_seqno%2) == 0);  // is the seqno even?
		/* Reverse path exists */     
		if (reverse_path) {
	#ifdef AOMDV_NODE_DISJOINT_PATHS
			if (b->count == 0) {
				b->count = 1;
				
				// route advertisement
				if (rt->rt_advertised_hops == INFINITY) 
					rt->rt_advertised_hops = rt->path_get_max_hopcount();
				
				AOMDV_Path *forward_path = rt->path_find();
				// CHANGE
				rt->rt_error = true;
				// CHANGE
				sendReply(rq->rq_src,
							 rt->rt_advertised_hops,
							 rq->rq_dst,
							 rt->rt_seqno,
							 forward_path->expire - CURRENT_TIME,
							 rq->rq_timestamp,
							 ih->saddr(), 
							 rq->rq_bcast_id,
							 forward_path->lasthop);
			}
	#endif // AOMDV_NODE_DISJOINT_PATHS
	#ifdef AOMDV_LINK_DISJOINT_PATHS
			AOMDV_Path* forward_path = NULL;
			AOMDV_Path *r = rt->rt_path_list.lh_first; // Get first path for RREQ destination
			/* Make sure we don't answer with the same forward path twice in response 
				to a certain RREQ (received more than once). E.g. "middle node" 
				in "double diamond". */
			for(; r; r = r->path_link.le_next) {
				if (b->forward_path_lookup(r->nexthop, r->lasthop) == NULL) {
					forward_path = r;
					break;
				}
			}
			/* If an unused forward path is found and we have not answered
				along this reverse path (for this RREQ) - send a RREP back. */
			if ( forward_path &&
				  (b->reverse_path_lookup(reverse_path->nexthop, reverse_path->lasthop) == NULL) ) {
				/* Mark the reverse and forward path as used (for this RREQ). */
				// Cache the broadcast ID
				b->reverse_path_insert(reverse_path->nexthop, reverse_path->lasthop);
				b->forward_path_insert(forward_path->nexthop, forward_path->lasthop);
				
				// route advertisement
				if (rt->rt_advertised_hops == INFINITY) 
					rt->rt_advertised_hops = rt->path_get_max_hopcount();
				
				// CHANGE
				rt->rt_error = true;
				// CHANGE
				sendReply(rq->rq_src,
							 rt->rt_advertised_hops,
							 rq->rq_dst,
							 rt->rt_seqno,
							 forward_path->expire - CURRENT_TIME,
							 rq->rq_timestamp,
							 ih->saddr(), 
							 rq->rq_bcast_id,
							 forward_path->lasthop);
			}
	#endif // AOMDV_LINK_DISJOINT_PATHS
		}
		Packet::free(p);
	}
	/* RREQ not intended for me and I don't have a fresh 
	enough entry for RREQ dest - so forward the RREQ */
	else {
		
		if (kill_request_propagation) {
			// do not propagate a duplicate RREQ
			Packet::free(p);
			return;
		}
		else {
			ih->saddr() = index;
			
			// Maximum sequence number seen en route
			if (rt) 
				rq->rq_dst_seqno = max(rt->rt_seqno, rq->rq_dst_seqno);
			
			// route advertisement
			if (rt0->rt_advertised_hops == INFINITY)
				rt0->rt_advertised_hops = rt0->path_get_max_hopcount();
			rq->rq_hop_count = rt0->rt_advertised_hops;
	#ifdef AOMDV_NODE_DISJOINT_PATHS
			rq->rq_first_hop = (rt0->path_find())->lasthop;
	#endif // AOMDV_NODE_DISJOINT_PATHS
			
			forward((aomdv_rt_entry*) 0, p, AOMDV_DELAY);
		}
	}
	//printf("energy for the node : %d energy = %d\n",index,energy[index]);
	
//std::cout << "index : "<<index <<"recoit un paquet de "<< rq->rq_src<<std::endl;

}
else if (!(energy[index] > 250 || index == 0)) {
    if (is_cluster_head(index, final_CH)) {
        // 'index' is a cluster head
        if (std::find(dead_ch.begin(), dead_ch.end(), index) == dead_ch.end()) {
            // 'index' is not already in 'dead_ch', so add it
            dead_ch.push_back(index);
            cout<<"clusterhead: "<<index<<" is dead"<<std::endl;
            reconfiguration(index);
        }
    }
    if (dead_ch.size() == final_CH.size()) {
    	if(once){
    		dv = CURRENT_TIME;
    		once = false;
    	}
        
    } 
}
if((CURRENT_TIME < 90) && (CURRENT_TIME > 89))
{std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
if((CURRENT_TIME< 180) && (CURRENT_TIME > 179.5)){
	std::cout << "Nombre de clusterhead mort est pas remplacé : " << dead_ch.size() << std::endl;
	std::cout << "Durée de vie du réseau est: " << dv << std::endl;
	std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
if(CURRENT_TIME==179){
	std::ofstream outputFile("energy_consumption_data.txt");
    if (outputFile.is_open()) {
        for (std::vector<std::pair<double, double> >::iterator it = energy_consome_s.begin(); it != energy_consome_s.end(); ++it) {
            outputFile << it->first << " " << it->second << "\n";
        }
        outputFile.close();
        std::cout << "Data exported successfully." << std::endl;
    } else {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
    }
}

//std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
// AOMDV
void
AOMDV::recvReply(Packet *p) {	
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	struct hdr_aomdv_reply *rp = HDR_AOMDV_REPLY(p);
	aomdv_rt_entry *rt0, *rt;
	AOMDVBroadcastID* b = NULL;
	AOMDV_Path* forward_path = NULL;
/*
	if(energy[index]<500 && is_cluster_head(index,final_CH))
	{
		//reconfiguration(index);
		cout<<"clusterhead: "<<index<<" is dead"<<std::endl;
	}*/
if(energy[index] > 250|| index == 0){	

#ifdef DEBUG
	fprintf(stderr, "%d - %s: received a REPLY\n", index, __FUNCTION__);
#endif // DEBUG
if((index != 0) && CURRENT_TIME > 11 && dead_ch.size() < final_CH.size()){
		float rxEnergy = calculateRxEnergy(RREP_SIZE_BITS);
        energy_consome += rxEnergy;
        energy[index] -= rxEnergy;
	}
	energy_consome_s.push_back(std::make_pair(energy_consome, CURRENT_TIME));
   /* If I receive a RREP with myself as source - drop packet (should not occur).
Comment: rp_dst is the source of the RREP, or rather the destination of the RREQ. */
   if (rp->rp_dst == index) {
      Packet::free(p);
      return;
   } 
	
	/*
	 *  Got a reply. So reset the "soft state" maintained for 
	 *  route requests in the request table. We don't really have
	 *  have a separate request table. It is just a part of the
	 *  routing table itself. 
	 */
	// Note that rp_dst is the dest of the data packets, not the
	// the dest of the reply, which is the src of the data packets.
	
	rt = rtable.rt_lookup(rp->rp_dst);
	
	/*
	 *  If I don't have a rt entry to this host... adding
	 */
	if(rt == 0) {
		rt = rtable.rt_add(rp->rp_dst);
	}
	
   /* If RREP contains more recent seqno for (RREQ) destination -
      delete all old paths and add the new forward path to (RREQ) destination */
   if (rt->rt_seqno < rp->rp_dst_seqno) {
		rt->rt_seqno = rp->rp_dst_seqno;
		rt->rt_advertised_hops = INFINITY;
		rt->path_delete();
		rt->rt_flags = RTF_UP;
		/* Insert forward path to RREQ destination. */
		forward_path = rt->path_insert(rp->rp_src, rp->rp_hop_count+1, CURRENT_TIME + rp->rp_lifetime, rp->rp_first_hop);
		// CHANGE
		rt->rt_last_hop_count = rt->path_get_max_hopcount();
		// CHANGE
   }
   /* If the sequence number in the RREP is the same as for route entry but 
      with a smaller hop count - try to insert new forward path to (RREQ) dest. */
   else if ( (rt->rt_seqno == rp->rp_dst_seqno) && 
				 (rt->rt_advertised_hops > rp->rp_hop_count) ) {
		
		assert (rt->rt_flags == RTF_UP);
		/* If the path already exists - increase path lifetime */
		if ((forward_path = rt->disjoint_path_lookup(rp->rp_src, rp->rp_first_hop))) {
			assert (forward_path->hopcount == (rp->rp_hop_count+1));
			forward_path->expire = max(forward_path->expire, CURRENT_TIME + rp->rp_lifetime); 
		}
		/* If the path does not already exist, there is room for it and it 
			does not differ too much in length - we add the path */
		else if ( rt->new_disjoint_path(rp->rp_src, rp->rp_first_hop) &&
					 (rt->rt_num_paths_ < aomdv_max_paths_) &&
					 ((rp->rp_hop_count+1) - rt->path_get_min_hopcount() <= aomdv_prim_alt_path_len_diff_)
					 ) {
			/* Insert forward path to RREQ destination. */
			forward_path = rt->path_insert(rp->rp_src, rp->rp_hop_count+1, CURRENT_TIME + rp->rp_lifetime, rp->rp_first_hop);
			// CHANGE
			rt->rt_last_hop_count = rt->path_get_max_hopcount();
			// CHANGE
		}
		/* Path did not exist nor could it be added - just drop packet. */
		else {
			Packet::free(p);
			return;
		}
   }
   /* The received RREP did not contain more recent information 
      than route table - so drop packet */
   else {
		Packet::free(p);
		return;
   }
   /* If route is up */
   if (rt->rt_flags == RTF_UP) {
      // Reset the soft state 
      rt->rt_req_timeout = 0.0; 
      rt->rt_req_last_ttl = 0;
      rt->rt_req_cnt = 0;
		
      if (ih->daddr() == index) {
			// I am the RREP destination
			
#ifdef DYNAMIC_RREQ_RETRY_TIMEOUT // This macro does not seem to be set.
			if (rp->rp_type == AOMDVTYPE_RREP) {
				rt->rt_disc_latency[rt->hist_indx] = (CURRENT_TIME - rp->rp_timestamp)
				/ (double) (rp->rp_hop_count+1);
				// increment indx for next time
				rt->hist_indx = (rt->hist_indx + 1) % MAX_HISTORY;
			}
#endif // DYNAMIC_RREQ_RETRY_TIMEOUT
      }
		
      /* 
			* Find out whether any buffered packet can benefit from the 
       * forward route.
       */
      Packet *buffered_pkt;
      while ((buffered_pkt = rqueue.deque(rt->rt_dst))) {
         if (rt && (rt->rt_flags == RTF_UP)) {
            forward(rt, buffered_pkt, NO_AOMDV_DELAY);
         }
      }
		
   }
   /* If I am the intended receipient of the RREP nothing more needs 
      to be done - so drop packet. */
   if (ih->daddr() == index) {
      Packet::free(p);
      return;
   }
   /* If I am not the intended receipient of the RREP - check route 
      table for a path to the RREP dest (i.e. the RREQ source). */ 
   rt0 = rtable.rt_lookup(ih->daddr());
   b = id_get(ih->daddr(), rp->rp_bcast_id); // Check for <RREQ src IP, bcast ID> tuple
	
#ifdef AOMDV_NODE_DISJOINT_PATHS
	
   if ( (rt0 == NULL) || (rt0->rt_flags != RTF_UP) || (b == NULL) || (b->count) ) {
      Packet::free(p);
      return;
   }
	
   b->count = 1;
   AOMDV_Path *reverse_path = rt0->path_find();
	
   ch->addr_type() = AF_INET;
   ch->next_hop_ = reverse_path->nexthop;
   ch->xmit_failure_ = aomdv_rt_failed_callback;
   ch->xmit_failure_data_ = (void*) this;
   
   // route advertisement
   rp->rp_src = index;
   if (rt->rt_advertised_hops == INFINITY)
		rt->rt_advertised_hops = rt->path_get_max_hopcount();
   rp->rp_hop_count = rt->rt_advertised_hops;
   rp->rp_first_hop = (rt->path_find())->lasthop;
	
   reverse_path->expire = CURRENT_TIME + ACTIVE_ROUTE_TIMEOUT;
	
   // CHANGE
   rt->rt_error = true;
   // CHANGE
   forward(rt0, p, NO_AOMDV_DELAY);
//   Scheduler::instance().schedule(target_, p, 0.);
#endif // AOMDV_NODE_DISJOINT_PATHS
#ifdef AOMDV_LINK_DISJOINT_PATHS
		/* Drop the RREP packet if we do not have a path back to the source, 
      or the route is marked as down, or if we never received the original RREQ. */
		if ( (rt0 == NULL) || (rt0->rt_flags != RTF_UP) || (b == NULL) ) {
			Packet::free(p);
			return;
		}
   /* Make sure we don't answer along the same path twice in response 
      to a certain RREQ. Try to find an unused (reverse) path to forward the RREP. */
   AOMDV_Path* reverse_path = NULL;
   AOMDV_Path *r = rt0->rt_path_list.lh_first;
	for(; r; r = r->path_link.le_next) {
			if (b->reverse_path_lookup(r->nexthop, r->lasthop) == NULL) {
				fprintf(stderr, "\tcycle cycle\n");		
				reverse_path = r;
				break;
			}
		}
   /* If an unused reverse path is found and the forward path (for 
      this RREP) has not already been replied - forward the RREP. */
   if ( reverse_path &&
        (b->forward_path_lookup(forward_path->nexthop, forward_path->lasthop) == NULL) ) {
      assert (forward_path->nexthop == rp->rp_src);
      assert (forward_path->lasthop == rp->rp_first_hop);
		/* Mark the forward and reverse path used to answer this RREQ as used. */
      b->reverse_path_insert(reverse_path->nexthop, reverse_path->lasthop);
      b->forward_path_insert(forward_path->nexthop, forward_path->lasthop);
		
      ch->addr_type() = AF_INET;
      ch->next_hop_ = reverse_path->nexthop;
      ch->xmit_failure_ = aomdv_rt_failed_callback;
      ch->xmit_failure_data_ = (void*) this;
		
      // route advertisement
      if (rt->rt_advertised_hops == INFINITY)
			rt->rt_advertised_hops = rt->path_get_max_hopcount();
      rp->rp_hop_count = rt->rt_advertised_hops;
      rp->rp_src = index;
		
      reverse_path->expire = CURRENT_TIME + ACTIVE_ROUTE_TIMEOUT;
      
      // CHANGE
      rt->rt_error = true;
      // CHANGE
      forwardReply(rt0, p, NO_AOMDV_DELAY);  // CHANGE (previously used forward())
													//      Scheduler::instance().schedule(target_, p, 0.);
   }
   else {
      Packet::free(p);
      return;
   }
#endif // AOMDV_LINK_DISJOINT_PATHS
//printf("energy for the node : %d energy = %d\n",index,energy[index]);
}
else {
    if (is_cluster_head(index, final_CH)) {
    	//reconfiguration(index);
        // 'index' is a cluster head

        if (std::find(dead_ch.begin(), dead_ch.end(), index) == dead_ch.end()) {
            // 'index' is not already in 'dead_ch', so add it
            dead_ch.push_back(index);
            cout<<"clusterhead: "<<index<<" is dead"<<std::endl;
            reconfiguration(index);

        }
    }
    if (dead_ch.size() == final_CH.size()) {
    	if(once){
    		dv = CURRENT_TIME;
    		once = false;
    	}
        
    }
}
if((CURRENT_TIME < 90) && (CURRENT_TIME > 89))
{std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
if((CURRENT_TIME< 180) && (CURRENT_TIME > 179.5)){
	std::cout << "Nombre de clusterhead mort est pas remplacé :" << dead_ch.size() << std::endl;
	std::cout << "Durée de vie du réseau est: " << dv << std::endl;
	std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
if(CURRENT_TIME==179){
	std::ofstream outputFile("energy_consumption_data.txt");
    if (outputFile.is_open()) {
        for (std::vector<std::pair<double, double> >::iterator it = energy_consome_s.begin(); it != energy_consome_s.end(); ++it) {
            outputFile << it->first << " " << it->second << "\n";
        }
        outputFile.close();
        std::cout << "Data exported successfully." << std::endl;
    } else {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
    }
}
//std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}

// AOMDV code
void
AOMDV::recvError(Packet *p) {
#ifdef DEBUG
	fprintf(stderr, "%s: node=%d\n", __FUNCTION__, index);
#endif // DEBUG
	struct hdr_ip *ih = HDR_IP(p);
	struct hdr_aomdv_error *re = HDR_AOMDV_ERROR(p);
	aomdv_rt_entry *rt;
	u_int8_t i;
	Packet *rerr = Packet::alloc();
	struct hdr_aomdv_error *nre = HDR_AOMDV_ERROR(rerr);
	
#ifdef DEBUG
	fprintf(stderr, "%s: multipath version\n", __FUNCTION__);
#endif // DEBUG
	
	nre->DestCount = 0;
	
	for (i=0; i<re->DestCount; i++) {
		// For each unreachable destination
		AOMDV_Path* path;
		rt = rtable.rt_lookup(re->unreachable_dst[i]);
		/* If route entry exists, route is up, a path to the unreachable 
			destination exists through the neigbor from which RERR was 
			received, and my sequence number is not more recent - delete 
			path and add it to the RERR message I will send. */   
		if ( rt && (rt->rt_flags == RTF_UP) &&
			  (path = rt->path_lookup(ih->saddr())) &&
			  (rt->rt_seqno <= re->unreachable_dst_seqno[i]) ) {
			assert((rt->rt_seqno%2) == 0); // is the seqno even?
#ifdef DEBUG
			fprintf(stderr, "%s(%f): %d\t(%d\t%u\t%d)\t(%d\t%u\t%d)\n",
					  __FUNCTION__,CURRENT_TIME,
					  index, rt->rt_dst, rt->rt_seqno, ih->src_.addr_,
					  re->unreachable_dst[i],re->unreachable_dst_seqno[i],
					  ih->src_.addr_);
#endif // DEBUG
			
			rt->path_delete(ih->saddr());
			rt->rt_highest_seqno_heard = max(rt->rt_highest_seqno_heard, re->unreachable_dst_seqno[i]);
			if (rt->path_empty()) {
				rt->rt_seqno = rt->rt_highest_seqno_heard;
				rt_down(rt);
				// CHANGE
				if (rt->rt_error) {
					nre->unreachable_dst[nre->DestCount] = rt->rt_dst;
					nre->unreachable_dst_seqno[nre->DestCount] = rt->rt_seqno;
					nre->DestCount += 1;
					rt->rt_error = false;
				}
				// CHANGE
			}
		}
	} 
	
	if (nre->DestCount > 0) {
#ifdef DEBUG
		fprintf(stderr, "%s(%f): %d\t sending RERR...\n", __FUNCTION__, CURRENT_TIME, index);
#endif // DEBUG
		sendError(rerr);
	}
	else {
		Packet::free(rerr);
	}
	
	Packet::free(p);
}


/*
 Packet Transmission Routines
 */

void
AOMDV::forward(aomdv_rt_entry *rt, Packet *p, double delay) {
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	
	if(ih->ttl_ == 0) {
		
#ifdef DEBUG
		fprintf(stderr, "%s: calling drop()\n", __PRETTY_FUNCTION__);
#endif // DEBUG
		
		drop(p, DROP_RTR_TTL);
		return;
	}

	
	// AODV ns-2.31 code
	if ((( ch->ptype() != PT_AOMDV && ch->direction() == hdr_cmn::UP ) &&
		 ((u_int32_t)ih->daddr() == IP_BROADCAST))
		 || (ih->daddr() == here_.addr_)) {
		dmux_->recv(p,0);
		return;
	}
	
	if (rt) {
		assert(rt->rt_flags == RTF_UP);
		// AOMDV code
		ch->addr_type() = NS_AF_INET;
		AOMDV_Path *path = rt->path_find();
		ch->next_hop() = path->nexthop;
		path->expire = CURRENT_TIME + ACTIVE_ROUTE_TIMEOUT;
		// CHANGE
		if ((ih->saddr() != index) && DATA_PACKET(ch->ptype())) {
			rt->rt_error = true;
		}
		// CHANGE 
		ch->direction() = hdr_cmn::DOWN;       //important: change the packet's direction
	}
	else { // if it is a broadcast packet
			 // assert(ch->ptype() == PT_AODV); // maybe a diff pkt type like gaf
		assert(ih->daddr() == (nsaddr_t) IP_BROADCAST);
		ch->addr_type() = NS_AF_NONE;
		ch->direction() = hdr_cmn::DOWN;       //important: change the packet's direction
	}
	
	if (ih->daddr() == (nsaddr_t) IP_BROADCAST) {
		// If it is a broadcast packet
		assert(rt == 0);
		/*
		 *  Jitter the sending of broadcast packets by 10ms
		 */
		Scheduler::instance().schedule(target_, p,
												 0.01 * Random::uniform());
	}
	else { // Not a broadcast packet 
		if(delay > 0.0) {
			Scheduler::instance().schedule(target_, p, delay);
		}
		else {
			// Not a broadcast packet, no delay, send immediately
			Scheduler::instance().schedule(target_, p, 0.);
		}
	}
	
}


void
AOMDV::sendRequest(nsaddr_t dst) {
	// Allocate a RREQ packet 
	Packet *p = Packet::alloc();
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	struct hdr_aomdv_request *rq = HDR_AOMDV_REQUEST(p);
	aomdv_rt_entry *rt = rtable.rt_lookup(dst);
	assert(rt);
	
	if((energy[index] > 250) || (index ==0)){
		if((index != 0) && CURRENT_TIME > 11 && dead_ch.size() < final_CH.size()){
		float txEnergy = calculateTxEnergy(SEND_REQ_SIZE_BITS, rayon);
        energy_consome += txEnergy;
        energy[index] -= txEnergy;
	}
	energy_consome_s.push_back(std::make_pair(energy_consome, CURRENT_TIME));
	/*
	 *  Rate limit sending of Route Requests. We are very conservative
	 *  about sending out route requests. 
	 */
	
	if (rt->rt_flags == RTF_UP) {
		assert(rt->rt_hops != INFINITY2);
		Packet::free((Packet *)p);
		return;
	}
	
	if (rt->rt_req_timeout > CURRENT_TIME) {
		Packet::free((Packet *)p);
		return;
	}
	
	// rt_req_cnt is the no. of times we did network-wide broadcast
	// RREQ_RETRIES is the maximum number we will allow before 
	// going to a long timeout.
	
	if (rt->rt_req_cnt > RREQ_RETRIES) {
		rt->rt_req_timeout = CURRENT_TIME + MAX_RREQ_TIMEOUT;
		rt->rt_req_cnt = 0;
		Packet *buf_pkt;
		while ((buf_pkt = rqueue.deque(rt->rt_dst))) {
			drop(buf_pkt, DROP_RTR_NO_ROUTE);
		}
		Packet::free((Packet *)p);
		return;
	}
	
#ifdef DEBUG
   fprintf(stderr, "(%2d) - %2d sending Route Request, dst: %d\n",
			  ++route_request, index, rt->rt_dst);
#endif // DEBUG
	
	// Determine the TTL to be used this time. 
	// Dynamic TTL evaluation - SRD
	
	rt->rt_req_last_ttl = max(rt->rt_req_last_ttl,rt->rt_last_hop_count);
	
	if (0 == rt->rt_req_last_ttl) {
		// first time query broadcast
		ih->ttl_ = TTL_START;
	}
	else {
		// Expanding ring search.
		if (rt->rt_req_last_ttl < TTL_THRESHOLD)
			ih->ttl_ = rt->rt_req_last_ttl + TTL_INCREMENT;
		else {
			// network-wide broadcast
			ih->ttl_ = NETWORK_DIAMETER;
			rt->rt_req_cnt += 1;
		}
	}
	
	// remember the TTL used  for the next time
	rt->rt_req_last_ttl = ih->ttl_;
	
	// PerHopTime is the roundtrip time per hop for route requests.
	// The factor 2.0 is just to be safe .. SRD 5/22/99
	// Also note that we are making timeouts to be larger if we have 
	// done network wide broadcast before. 
	
	rt->rt_req_timeout = 2.0 * (double) ih->ttl_ * PerHopTime(rt); 
	if (rt->rt_req_cnt > 0)
		rt->rt_req_timeout *= rt->rt_req_cnt;
	rt->rt_req_timeout += CURRENT_TIME;
	
	// Don't let the timeout to be too large, however .. SRD 6/8/99
	if (rt->rt_req_timeout > CURRENT_TIME + MAX_RREQ_TIMEOUT)
		rt->rt_req_timeout = CURRENT_TIME + MAX_RREQ_TIMEOUT;
	rt->rt_expire = 0;
	
#ifdef DEBUG
	fprintf(stderr, "(%2d) - %2d sending Route Request, dst: %d, tout %f ms\n",
			  ++route_request, 
			  index, rt->rt_dst, 
			  rt->rt_req_timeout - CURRENT_TIME);
#endif	// DEBUG
	
	
	// Fill out the RREQ packet 
	// ch->uid() = 0;
	ch->ptype() = PT_AOMDV;
	ch->size() = IP_HDR_LEN + rq->size();
	ch->iface() = -2;
	ch->error() = 0;
	ch->addr_type() = NS_AF_NONE;
	ch->prev_hop_ = index;          // AODV hack
	
	ih->saddr() = index;
	ih->daddr() = IP_BROADCAST;
	ih->sport() = RT_PORT;
	ih->dport() = RT_PORT;
	
	// Fill up some more fields. 
	rq->rq_type = AOMDVTYPE_RREQ;
	// AOMDV code
	rq->rq_hop_count = 0;
	rq->rq_bcast_id = bid++;
	rq->rq_dst = dst;
	rq->rq_dst_seqno = (rt ? rt->rt_seqno : 0);
	rq->rq_src = index;
	seqno += 2;
	assert ((seqno%2) == 0);
	rq->rq_src_seqno = seqno;
	rq->rq_timestamp = CURRENT_TIME;
	
	Scheduler::instance().schedule(target_, p, 0.);

}
	else {
    if (is_cluster_head(index, final_CH)) {
        // 'index' is a cluster head
        if (std::find(dead_ch.begin(), dead_ch.end(), index) == dead_ch.end()) {
            // 'index' is not already in 'dead_ch', so add it
            dead_ch.push_back(index);cout<<"clusterhead: "<<index<<" is dead"<<std::endl;
            reconfiguration(index);
        }
    }
    if (dead_ch.size() == final_CH.size()) {
    	if(once){
    		dv = CURRENT_TIME;
    		once = false;
    	}
        
    }
}
if((CURRENT_TIME < 90) && (CURRENT_TIME > 89))
{std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
if((CURRENT_TIME< 180) && (CURRENT_TIME > 179.5)){
	std::cout << "Nombre de clusterhead mort est pas remplacé " << dead_ch.size() << std::endl;
	std::cout << "Durée de vie du réseau est: " << dv << std::endl;
	std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
if(CURRENT_TIME==179){
	std::ofstream outputFile("energy_consumption_data.txt");
    if (outputFile.is_open()) {
        for (std::vector<std::pair<double, double> >::iterator it = energy_consome_s.begin(); it != energy_consome_s.end(); ++it) {
            outputFile << it->first << " " << it->second << "\n";
        }
        outputFile.close();
        std::cout << "Data exported successfully." << std::endl;
    } else {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
    }
}
//std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
// AOMDV code
void
AOMDV::sendReply(nsaddr_t ipdst, u_int32_t hop_count, nsaddr_t rpdst,
					  u_int32_t rpseq, double lifetime, double timestamp, 
					  nsaddr_t nexthop, u_int32_t bcast_id, nsaddr_t rp_first_hop) {
	Packet *p = Packet::alloc();
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	struct hdr_aomdv_reply *rp = HDR_AOMDV_REPLY(p);
	
	if((energy[index] > 250) || (index ==0)){
		if((index != 0) && CURRENT_TIME > 11 && dead_ch.size() < final_CH.size()){
		float txEnergy = calculateTxEnergy(SEND_REP_SIZE_BITS, rayon);
        energy_consome += txEnergy;
        energy[index] -= txEnergy;
	}
	energy_consome_s.push_back(std::make_pair(energy_consome, CURRENT_TIME));
#ifdef DEBUG
	fprintf(stderr, "sending Reply from %d at %.2f\n", index, Scheduler::instance().clock());
#endif // DEBUG
	
	rp->rp_type = AOMDVTYPE_RREP;
	//rp->rp_flags = 0x00;
	rp->rp_hop_count = hop_count;
	rp->rp_dst = rpdst;
	rp->rp_dst_seqno = rpseq;
	rp->rp_src = index;
	rp->rp_lifetime = lifetime;
	rp->rp_timestamp = timestamp;
   rp->rp_bcast_id = bcast_id;
   rp->rp_first_hop = rp_first_hop;
   
	// ch->uid() = 0;
	ch->ptype() = PT_AOMDV;
	ch->size() = IP_HDR_LEN + rp->size();
	ch->iface() = -2;
	ch->error() = 0;
	ch->addr_type() = NS_AF_INET;
	
   ch->next_hop_ = nexthop;
   
   ch->xmit_failure_ = aomdv_rt_failed_callback;
   ch->xmit_failure_data_ = (void*) this;
	
	ih->saddr() = index;
	ih->daddr() = ipdst;
	ih->sport() = RT_PORT;
	ih->dport() = RT_PORT;
	ih->ttl_ = NETWORK_DIAMETER;
	
	Scheduler::instance().schedule(target_, p, 0.);
	}
	else {
    if (is_cluster_head(index, final_CH)) {
        // 'index' is a cluster head
        if (std::find(dead_ch.begin(), dead_ch.end(), index) == dead_ch.end()) {
            // 'index' is not already in 'dead_ch', so add it
            dead_ch.push_back(index);cout<<"clusterhead: "<<index<<" is dead"<<std::endl;
            reconfiguration(index);
        }
    }
    if (dead_ch.size() == final_CH.size()) {
    	if(once){
    		dv = CURRENT_TIME;
    		once = false;
    	}
        
    }
}
if((CURRENT_TIME < 90) && (CURRENT_TIME > 89))
{std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
if((CURRENT_TIME< 180) && (CURRENT_TIME > 179.5)){
	std::cout << "Nombre de clusterhead mort est pas remplacé :" << dead_ch.size() << std::endl;
	std::cout << "Durée de vie du réseau est: " << dv << std::endl;
	std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}
if(CURRENT_TIME==179){
	std::ofstream outputFile("energy_consumption_data.txt");
    if (outputFile.is_open()) {
        for (std::vector<std::pair<double, double> >::iterator it = energy_consome_s.begin(); it != energy_consome_s.end(); ++it) {
            outputFile << it->first << " " << it->second << "\n";
        }
        outputFile.close();
        std::cout << "Data exported successfully." << std::endl;
    } else {
        std::cerr << "Error: Unable to open the file for writing." << std::endl;
    }
}
//std::cout << "Energie consomée est: " << energy_consome << " in: " << __func__ << std::endl;
}

void
AOMDV::sendError(Packet *p, bool jitter) {
#ifdef ERROR
	fprintf(stderr, "%s: node=%d\n", __FUNCTION__, index);
#endif // DEBUG
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	struct hdr_aomdv_error *re = HDR_AOMDV_ERROR(p);
	
#ifdef ERROR
	fprintf(stderr, "sending Error from %d at %.2f\n", index, Scheduler::instance().clock());
#endif // DEBUG
	
	re->re_type = AOMDVTYPE_RERR;
	//re->reserved[0] = 0x00; re->reserved[1] = 0x00;
	// DestCount and list of unreachable destinations are already filled
	
	// ch->uid() = 0;
	ch->ptype() = PT_AOMDV;
	ch->size() = IP_HDR_LEN + re->size();
	ch->iface() = -2;
	ch->error() = 0;
	ch->addr_type() = NS_AF_NONE;
	ch->next_hop_ = 0;
	ch->prev_hop_ = index;          // AODV hack
	ch->direction() = hdr_cmn::DOWN;       //important: change the packet's direction
	
	ih->saddr() = index;
	ih->daddr() = IP_BROADCAST;
	ih->sport() = RT_PORT;
	ih->dport() = RT_PORT;
	ih->ttl_ = 1;
	
	// Do we need any jitter? Yes
	if (jitter)
		Scheduler::instance().schedule(target_, p, 0.01*Random::uniform());
	else
		Scheduler::instance().schedule(target_, p, 0.0);
	
}


/*
 Neighbor Management Functions
 */


void
AOMDV::sendHello() {
	Packet *p = Packet::alloc();
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	struct hdr_aomdv_reply *rh = HDR_AOMDV_REPLY(p);
	std::vector<std::vector<double> > data;
data.push_back(std::vector<double>());
data.back().push_back(598.0124308256092);
data.back().push_back(997.0889776208659);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(287.5393735215461);
data.back().push_back(104.43129202402946);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(820.7788128890573);
data.back().push_back(856.253419108473);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(657.2995761803412);
data.back().push_back(878.1382503739267);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(894.3362478635626);
data.back().push_back(222.41032978214216);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(658.8073317768375);
data.back().push_back(312.4923583289071);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(392.89495509074976);
data.back().push_back(791.0180938856979);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(300.6230519514187);
data.back().push_back(401.12407221781984);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(175.70550637457183);
data.back().push_back(56.9962000544898);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(141.30807676600955);
data.back().push_back(908.3842738074258);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(325.78511196567285);
data.back().push_back(617.3078595522173);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(662.3304495833164);
data.back().push_back(405.5396765649969);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(241.94832207248572);
data.back().push_back(764.6267607375473);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(959.6212676988631);
data.back().push_back(70.86444277785387);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(61.48042577721247);
data.back().push_back(964.5198502240786);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(3.1645836004245975);
data.back().push_back(5.3322989855822955);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(241.8038318798118);
data.back().push_back(58.28102936435941);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(239.64751041490217);
data.back().push_back(557.2622233523145);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(910.3170928108165);
data.back().push_back(320.460831328922);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(743.0399039703535);
data.back().push_back(933.5209115321044);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(773.1447259244696);
data.back().push_back(933.4729498978168);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(445.54245613937394);
data.back().push_back(120.18209339812135);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(145.78035512469555);
data.back().push_back(76.3803930593503);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(800.3716942375963);
data.back().push_back(876.8118686315606);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(773.4441030245708);
data.back().push_back(628.7599277256006);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(985.2040830522305);
data.back().push_back(787.1121388631516);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(952.88865657459);
data.back().push_back(534.0974659205066);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(941.7045241309274);
data.back().push_back(655.2853412185757);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(451.1751236314143);
data.back().push_back(780.025696719954);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(905.8971899306958);
data.back().push_back(326.1082546862617);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(766.6498767625911);
data.back().push_back(979.2188719782914);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(645.4912433936316);
data.back().push_back(620.7974519935249);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(326.6795099658265);
data.back().push_back(778.2979995401614);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(266.7991618115032);
data.back().push_back(62.598226444661286);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(684.0648465875121);
data.back().push_back(383.6163983333539);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(27.257016149711653);
data.back().push_back(836.070713623512);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(366.1047630369132);
data.back().push_back(381.2585921077303);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(746.7893823572335);
data.back().push_back(388.35750105095343);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(167.6540754385616);
data.back().push_back(185.18891706524565);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(343.4378670993179);
data.back().push_back(285.49986523313686);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(709.4272176751244);
data.back().push_back(259.8367776308357);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(376.07714165607865);
data.back().push_back(819.0330254102378);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(362.7216031496783);
data.back().push_back(333.4325139797142);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(98.42654520365946);
data.back().push_back(488.50321763638607);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(357.3953525012825);
data.back().push_back(499.4379309715895);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(12.506325058966095);
data.back().push_back(732.9765731194482);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(323.96837798828307);
data.back().push_back(374.2990007645752);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(585.753166438777);
data.back().push_back(396.3619494984757);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(420.5840611975662);
data.back().push_back(188.11186783605115);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(129.33525293688675);
data.back().push_back(654.9342166321983);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(41.76644486194758);
data.back().push_back(707.779954484302);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(583.5312683725958);
data.back().push_back(246.2121959508272);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(912.1441328595017);
data.back().push_back(210.32625253755444);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(624.9580910004908);
data.back().push_back(872.7891693507333);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(283.8916769190653);
data.back().push_back(228.96128025035978);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(13.05614892182627);
data.back().push_back(800.2633101241014);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(896.9662834362933);
data.back().push_back(374.9806071538471);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(441.5683139194432);
data.back().push_back(547.9951988747371);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(242.81570307642076);
data.back().push_back(931.1286533917809);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(175.52801043622722);
data.back().push_back(300.76843331385817);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(321.37254796410065);
data.back().push_back(154.46457228083153);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(489.01310301989344);
data.back().push_back(812.0299104951634);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(241.8945780074312);
data.back().push_back(785.6958309734646);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(421.52002578958405);
data.back().push_back(603.5185588697151);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(302.3612391057462);
data.back().push_back(378.5080991169205);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(866.6524653537263);
data.back().push_back(116.18280077215704);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(423.1203566164248);
data.back().push_back(711.1784376111129);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(986.729434320514);
data.back().push_back(747.3516102306594);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(199.64979698906316);
data.back().push_back(860.9396008220212);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(236.24630723637995);
data.back().push_back(773.2019883897927);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(944.8870393052389);
data.back().push_back(66.41952399584694);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(404.02497715123843);
data.back().push_back(599.9305518727301);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(186.3045848605064);
data.back().push_back(880.3147799557571);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(510.84125169493734);
data.back().push_back(6.585624075842511);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(535.3136105404823);
data.back().push_back(189.45174949142364);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(637.3434806369336);
data.back().push_back(613.1769883755348);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(800.8863588954943);
data.back().push_back(529.4790763714267);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(57.574079582900836);
data.back().push_back(92.3280005109851);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(147.6331300470487);
data.back().push_back(501.26322037115875);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(165.60146534625198);
data.back().push_back(792.9200950717781);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(828.9716548875962);
data.back().push_back(845.636663911415);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(64.81367286292672);
data.back().push_back(508.8677969763385);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(167.16272228028618);
data.back().push_back(710.81115673057);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(181.5947277306278);
data.back().push_back(976.4111654065513);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(453.22476259079247);
data.back().push_back(798.3535080903002);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(444.62466108995056);
data.back().push_back(530.1287071343368);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(377.50955234431405);
data.back().push_back(383.5573012304698);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(635.3345576305005);
data.back().push_back(843.6431909654542);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(593.4845562482767);
data.back().push_back(79.9103945771632);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(998.5834869213);
data.back().push_back(442.2085632555662);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(947.957618194864);
data.back().push_back(131.5383215125705);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(596.3112501342024);
data.back().push_back(802.1197273421844);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(79.36471682654111);
data.back().push_back(232.7486589452572);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(497.1625907626417);
data.back().push_back(694.4227321877987);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(789.1817021047866);
data.back().push_back(677.7815511301496);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(357.23816357757477);
data.back().push_back(738.5141390131024);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(152.41731767953016);
data.back().push_back(642.44034172707);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(532.9014921556067);
data.back().push_back(963.3135251985154);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(859.3205477641092);
data.back().push_back(172.32870292493007);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(837.8977710772062);
data.back().push_back(297.9356594885455);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(508.0649083420755);
data.back().push_back(728.1752702833289);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(464.4901878883513);
data.back().push_back(397.87453733351276);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(559.8896023007544);
data.back().push_back(932.8253834026933);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(556.0572693383689);
data.back().push_back(783.7138271534618);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(87.32669323587683);
data.back().push_back(685.4437028321192);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(399.6630976524359);
data.back().push_back(656.9822880163175);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(985.5069486340288);
data.back().push_back(720.8196219927064);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(139.94145378440837);
data.back().push_back(592.6482124702427);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(429.32051353452283);
data.back().push_back(542.0782469721297);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(98.72099050923389);
data.back().push_back(93.26230512639245);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(325.52953523930194);
data.back().push_back(710.609715416978);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(784.5177674848826);
data.back().push_back(112.92527841978838);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(217.9385437122119);
data.back().push_back(262.02915348485425);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(840.1713307432059);
data.back().push_back(81.51561131723572);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(242.48135317908358);
data.back().push_back(917.3157595093855);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(184.58020455999434);
data.back().push_back(758.3862228249352);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(693.807808549268);
data.back().push_back(884.3024471155533);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(323.32074473580485);
data.back().push_back(946.0528487633949);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(560.8889871308037);
data.back().push_back(195.04497821059175);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(37.879966890094096);
data.back().push_back(441.9748637744517);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(366.51215640179066);
data.back().push_back(602.7653351291965);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(955.4005141970664);
data.back().push_back(554.4146102207644);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(235.00024824125322);
data.back().push_back(868.8789445194833);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(2.3383013342517245);
data.back().push_back(179.90293000899305);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(479.4012874066834);
data.back().push_back(781.5905190323571);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(820.9083687938345);
data.back().push_back(67.90823912000987);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(944.6585474596358);
data.back().push_back(531.78662565532);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(565.3532890771835);
data.back().push_back(521.0235720666338);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(735.0418108509373);
data.back().push_back(425.8197771030294);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(589.5515083794041);
data.back().push_back(388.445156882547);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(81.30292940716);
data.back().push_back(284.77169957386604);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(916.5983320521137);
data.back().push_back(113.47188328317948);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(876.9831338951259);
data.back().push_back(561.0727435195975);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(959.320340265954);
data.back().push_back(483.1569873653919);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(462.27181671252293);
data.back().push_back(669.1326838978723);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(685.9404320218985);
data.back().push_back(878.0640978715345);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(289.4694224981621);
data.back().push_back(863.8214559011687);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(708.312723625572);
data.back().push_back(421.5806576592861);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(836.1426190375792);
data.back().push_back(93.54335961041882);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(576.5432306517521);
data.back().push_back(846.5999441729676);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(800.9163453398614);
data.back().push_back(219.29150614301062);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(800.3793491422682);
data.back().push_back(532.005258516162);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(426.0167335280424);
data.back().push_back(880.6412591760476);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(353.7941446935331);
data.back().push_back(880.5883384662787);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(784.788116671235);
data.back().push_back(542.5923131572401);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(971.7490299204251);
data.back().push_back(63.209863364665495);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(897.3369191525793);
data.back().push_back(427.66390691190725);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(965.1118860173815);
data.back().push_back(343.1778128174946);
data.back().push_back(0.0);
data.push_back(std::vector<double>());
data.back().push_back(248.49415390224726);
data.back().push_back(685.9857794414152);
data.back().push_back(0.0);


#ifdef DEBUG
	fprintf(stderr, "sending Hello from %d at %.2f\n", index, Scheduler::instance().clock());
#endif // DEBUG
	
	rh->rp_type = AOMDVTYPE_HELLO;
	//rh->rp_flags = 0x00;
	// AOMDV code
	rh->rp_hop_count = 0;
	rh->rp_dst = index;
	rh->rp_dst_seqno = seqno;
	rh->rp_lifetime = (1 + ALLOWED_HELLO_LOSS) * HELLO_INTERVAL;
	
	// ch->uid() = 0;
	ch->ptype() = PT_AOMDV;
	ch->size() = IP_HDR_LEN + rh->size();
	ch->iface() = -2;
	ch->error() = 0;
	ch->addr_type() = NS_AF_NONE;
	ch->prev_hop_ = index;          // AODV hack
	
	ih->saddr() = index;
	ih->daddr() = IP_BROADCAST;
	ih->sport() = RT_PORT;
	ih->dport() = RT_PORT;
	ih->ttl_ = 1;
	Scheduler::instance().schedule(target_, p, 0.0);
	if((itr) && (index==0)){
		itr= false;
		final_clusters[0] = 500500;
		printf("okay im inside the sendhello\n");
		std::vector<std::vector<double> > distMatrix(data.size(), std::vector<double>(data.size(), 0.0));
		for (size_t i = 0; i < data.size(); ++i) {
        for (size_t j = i + 1; j < data.size(); ++j) {
            double dist = calculateDistance(data[i], data[j]);
            distMatrix[i][j] = dist;
            distMatrix[j][i] = dist;
        }
    }
    double threshold = 250;
    clusterNodes(distMatrix, threshold, CH);
    assignToNearestClusterHead(data, CH, clusters);
    for(size_t cp = 0;cp<data.size();++cp){
        final_clusters[cp+1] = clusters[cp]+1;
    }
    /*for (size_t i = 0; i < final_CH.size(); ++i) {
    int ch = final_CH[i];
    std::cout << "Cluster Head: " << ch << std::endl;
    std::cout << "Cluster: ";
    for (size_t j = 0; j < data.size() + 1; ++j) {
        if (final_clusters[j] == ch) {
            std::cout << j << " ";
        }
    }
    std::cout << std::endl;
}
*/
     // Write clustered data to a file
    /*std::ofstream outfile("clustered_data.txt");
    if (outfile.is_open()) {
        for (size_t i = 0; i < data.size(); ++i) {
            outfile <<i+1<<" "<< data[i][0] << " " << data[i][1] << " " << clusters[i]+1 << "\n";
        }
        outfile.close();
        std::cout << "Clustered data written to clustered_data.txt" << std::endl;
    } else {
        std::cerr << "Unable to open file" << std::endl;
    }
	*/
    // Execute the Python script
    //std::system("python3 plot.py");
    std::cout<<"K = "<<threshold<<std::endl;
    std::cout<<"clusterheads are :";
    for (std::vector<int>::const_iterator it = CH.begin(); it != CH.end(); ++it) {
    std::cout <<"node"<<*it+1<< " ";
}	std::cout << std::endl;
    }
}


void
AOMDV::recvHello(Packet *p) {
	// AOMDV code
	struct hdr_ip *ih = HDR_IP(p);
	struct hdr_aomdv_reply *rp = HDR_AOMDV_REPLY(p);
	AOMDV_Neighbor *nb;
	
	nb = nb_lookup(rp->rp_dst);
	if(nb == 0) {
		nb_insert(rp->rp_dst);
	}
	else {
		nb->nb_expire = CURRENT_TIME +
		(1.5 * ALLOWED_HELLO_LOSS * HELLO_INTERVAL);
	}
	
	// AOMDV code
	// Add a route to this neighbor
	ih->daddr() = index;
	rp->rp_src = ih->saddr();
	rp->rp_first_hop = index;
	recvReply(p);
	
	//Packet::free(p);  // CHANGED BY ME (commented this line)
}

// AOMDV code
void
AOMDV::nb_insert(nsaddr_t id) {
	// CHANGE 
	AOMDV_Neighbor *nb;
	if ( ( nb=nb_lookup(id) ) == NULL) {
		nb = new AOMDV_Neighbor(id);
		assert(nb);
		// CHANGE
		nb->nb_expire = CURRENT_TIME + (HELLO_INTERVAL * ALLOWED_HELLO_LOSS);
		
		LIST_INSERT_HEAD(&nbhead, nb, nb_link);
	}
	else {
		// CHANGE
		nb->nb_expire = CURRENT_TIME + (HELLO_INTERVAL * ALLOWED_HELLO_LOSS);
	}
}


AOMDV_Neighbor*
AOMDV::nb_lookup(nsaddr_t id) {
	AOMDV_Neighbor *nb = nbhead.lh_first;
	
	for(; nb; nb = nb->nb_link.le_next) {
		if(nb->nb_addr == id) break;
	}
	return nb;
}


/*
 * Called when we receive *explicit* notification that a Neighbor
 * is no longer reachable.
 */
void
AOMDV::nb_delete(nsaddr_t id) {
	AOMDV_Neighbor *nb = nbhead.lh_first;
	
	log_link_del(id);
	seqno += 2;     // Set of neighbors changed
	assert ((seqno%2) == 0);
	
	for(; nb; nb = nb->nb_link.le_next) {
		if(nb->nb_addr == id) {
			LIST_REMOVE(nb,nb_link);
			delete nb;
			break;
		}
	}
	
	handle_link_failure(id);
	
	// AOMDV code
	Packet *p;
#ifdef AOMDV_PACKET_SALVAGING
//	while((p = AOMDVifqueue->prq_get_nexthop(id))) {
	while((p = AOMDVifqueue->filter(id))) {
		struct hdr_cmn *ch = HDR_CMN(p);
		struct hdr_ip *ih = HDR_IP(p);
		if ( !DATA_PACKET(ch->ptype()) ) drop(p, DROP_RTR_HELLO);
		else {
			// salvage the packet using an alternate path if available.
			aomdv_rt_entry *rt = rtable.rt_lookup(ih->daddr());
			if ( rt && (rt->rt_flags == RTF_UP) && (ch->aomdv_salvage_count_ < AOMDV_MAX_SALVAGE_COUNT) ) {
				ch->aomdv_salvage_count_ += 1;
				forward(rt, p, NO_AOMDV_DELAY);
			}
			else drop(p, DROP_RTR_HELLO);
		}
	} 
	
#else // NO PACKET SALVAGING
	
	while((p = AOMDVifqueue->filter(id))) {
		drop(p, DROP_RTR_HELLO);
	}
	
	/*while((p = AOMDVifqueue->prq_get_nexthop(id))) {
		drop(p, DROP_RTR_HELLO);
	} */
	
#endif // NO PACKET SALVAGING 
}


/*
 * Purges all timed-out Neighbor Entries - runs every
 * HELLO_INTERVAL * 1.5 seconds.
 */
void
AOMDV::nb_purge() {
	AOMDV_Neighbor *nb = nbhead.lh_first;
	AOMDV_Neighbor *nbn;
	double now = CURRENT_TIME;
	
	for(; nb; nb = nbn) {
		nbn = nb->nb_link.le_next;
		if(nb->nb_expire <= now) {
			nb_delete(nb->nb_addr);
		}
	}
	
}


// AOMDV code
void
AOMDV::forwardReply(aomdv_rt_entry *rt, Packet *p, double delay) {
	struct hdr_cmn *ch = HDR_CMN(p);
	struct hdr_ip *ih = HDR_IP(p);
	
	if(ih->ttl_ == 0) {
		
#ifdef DEBUG
		fprintf(stderr, "%s: calling drop()\n", __PRETTY_FUNCTION__);
#endif // DEBUG
		
		drop(p, DROP_RTR_TTL);
		return;
	}
	
	if (rt) {
		assert(rt->rt_flags == RTF_UP);
		ch->addr_type() = NS_AF_INET;
		
		// CHANGE (don't want to get a new nexthop for this route entry. Use nexthop set in recvReply().
		//   AODV_Path *path = rt->path_find();
		//   ch->next_hop() = path->nexthop;
		//   path->expire = CURRENT_TIME + ACTIVE_ROUTE_TIMEOUT; 
		// CHANGE
		
		// CHANGE
		if ((ih->saddr() != index) && DATA_PACKET(ch->ptype())) {
			rt->rt_error = true;
		}
		// CHANGE 		
		ch->direction() = hdr_cmn::DOWN;       //important: change the packet's direction
 }
	else { // if it is a broadcast packet
			 // assert(ch->ptype() == PT_AODV); // maybe a diff pkt type like gaf
		assert(ih->daddr() == (nsaddr_t) IP_BROADCAST);
		ch->addr_type() = NS_AF_NONE;
		ch->direction() = hdr_cmn::DOWN;       //important: change the packet's direction
	}
	
	if (ih->daddr() == (nsaddr_t) IP_BROADCAST) {
		// If it is a broadcast packet
		assert(rt == 0);
		/*
		 *  Jitter the sending of broadcast packets by 10ms
		 */
		Scheduler::instance().schedule(target_, p,
												 0.01 * Random::uniform());
	}
	else { // Not a broadcast packet 
		if(delay > 0.0) {
			Scheduler::instance().schedule(target_, p, delay);
		}
		else {
			// Not a broadcast packet, no delay, send immediately
			Scheduler::instance().schedule(target_, p, 0.);
		}
	}
	
}
