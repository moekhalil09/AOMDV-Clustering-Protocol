#======================================================================
# Parametros de la linea de comandos
#======================================================================
global argv arg0
set opt(cp) [lindex $argv 2];#scenarios de traffic 
set opt(sc) [lindex $argv 3];# scenarios de mobilité 
set opt(tr) [lindex $argv 0];# fichier trace
set opt(na) [lindex $argv 1];# fichier nam
 
puts "----------------------------"
puts $opt(tr)
puts $opt(na)
puts $opt(cp)
puts $opt(sc)
puts "----------------------------"

# Define options
# ======================================================================
set opt(chan)           Channel/WirelessChannel  ;# channel type
set opt(prop)           Propagation/Nakagami 	;# radio-propagation model
set opt(netif)          Phy/WirelessPhyExt      ;# network interface type
set opt(mac)            Mac/802_11Ext           ;# MAC type

set opt(ifq)	Queue/DropTail/PriQueue
set opt(ll)		LL
set opt(ant)        Antenna/OmniAntenna

set opt(x)		2000.06		;# the width (x-axis) of the simulation
set opt(y)		1114.25		;# the width (x-axis) of the simulation
set opt(ifqlen)	50	      ;# max packet in ifq
set opt(seed)	0.0
set opt(tr)		aodmv.tr    ;# fichier contenant les resultat                        <-- tu dois changer ça 
set opt(adhocRouting)   AOMDV
set opt(nn)            150       ;# how many nodes are simulated
set opt(cp)		"scene/scene"   ;# scenario de transmission                   <-- change to path of cbr    
set opt(sc)		"cbr/scenecbr.tcl" # <-- change to path of scene
set opt(stop)		180.0		;# simulation time



#======================================================================
# Main Program
# ======================================================================
# Initialize Global Variables
# create simulator instance
set ns_		[new Simulator]

# set wireless channel, radio-model and topography objects
set wtopo	[new Topography]

# create trace object for ns and nam
set tracefd	[open $opt(tr) w]
$ns_ trace-all $tracefd
# use new trace file format
$ns_ use-newtrace 


 set namtrace [open aodv10.nam w]
      $ns_ namtrace-all-wireless $namtrace  $opt(x) $opt(y)
# ###########################################################
# define topology
$wtopo load_flatgrid $opt(x) $opt(y)




# ======================================================================
# Define IEEE 802.11P
# ======================================================================
Phy/WirelessPhyExt set CSThresh_                3.162e-12   ;#-85 dBm Wireless interface sensitivity (sensitivity defined in the standard)
Phy/WirelessPhyExt set Pt_                      0.005         
Phy/WirelessPhyExt set freq_                    5.85e+9; #5.9e+9
Phy/WirelessPhyExt set noise_floor_             1.26e-13    ;#-99 dBm for 10MHz bandwidth
Phy/WirelessPhyExt set L_                       1.0         ;#default radio circuit gain/loss
Phy/WirelessPhyExt set PowerMonitorThresh_      6.310e-14   ;#-102dBm power monitor  sensitivity
Phy/WirelessPhyExt set HeaderDuration_          0.000040    ;#40 us
Phy/WirelessPhyExt set BasicModulationScheme_   0
Phy/WirelessPhyExt set PreambleCaptureSwitch_   1
Phy/WirelessPhyExt set DataCaptureSwitch_       0
Phy/WirelessPhyExt set SINR_PreambleCapture_    2.5118;     ;# 4 dB
Phy/WirelessPhyExt set SINR_DataCapture_        100.0;      ;# 10 dB
Phy/WirelessPhyExt set trace_dist_              1e6         ;# PHY trace until distance of 1 Mio. km ("infinty")
Phy/WirelessPhyExt set PHY_DBG_                 0
Phy/WirelessPhyExt set bandwidth_ 70e6; #DEL PAPER

Mac/802_11Ext set CWMin_                        15
Mac/802_11Ext set CWMax_                        1023
Mac/802_11Ext set SlotTime_                     0.000013
Mac/802_11Ext set SIFS_                         0.000032
Mac/802_11Ext set ShortRetryLimit_              7
Mac/802_11Ext set LongRetryLimit_               4
Mac/802_11Ext set HeaderDuration_               0.000040
Mac/802_11Ext set SymbolDuration_               0.000008
Mac/802_11Ext set BasicModulationScheme_        0
Mac/802_11Ext set use_802_11a_flag_             true
Mac/802_11Ext set RTSThreshold_                 2346
Mac/802_11Ext set MAC_DBG                       0
# ======================================================================


# ======================================================================
# Urban propagation model: Nakagami
# ======================================================================
Propagation/Nakagami set  use_nakagami_dist_ true 
Propagation/Nakagami set gamma0_ 2.0 
Propagation/Nakagami set gamma1_ 2.0 
Propagation/Nakagami set gamma2_ 2.0 
Propagation/Nakagami set d0_gamma_ 200 
Propagation/Nakagami set d1_gamma_ 500 
Propagation/Nakagami set m0_  1.0 
Propagation/Nakagami set m1_  1.0 
Propagation/Nakagami set m2_  1.0 
Propagation/Nakagami set d0_m_ 80 
Propagation/Nakagami set d1_m_ 200 
# ======================================================================











# Create God
set god_ [create-god $opt(nn)]


# ======================================================================
# configure mobile nodes
# ======================================================================
# define how node should be created
#global node setting
$ns_ node-config -adhocRouting $opt(adhocRouting) \
			 -llType $opt(ll) \
		 -macType $opt(mac) \
		 -ifqType $opt(ifq) \
		 -ifqLen $opt(ifqlen) \
		 -antType $opt(ant) \
		 -propType $opt(prop) \
		 -phyType $opt(netif) \
		 -channelType $opt(chan) \
		 -topoInstance $wtopo \
		 -agentTrace ON \
                 -routerTrace ON \
                 -macTrace OFF 

				 
				 
#  Create the specified number of nodes [$opt(nn)] and "attach" them
#  to the channel. 
#for {set i 0} {$i < $opt(nn)-3 } {incr i} {
for {set i 0} {$i < $opt(nn) } {incr i} {
	set node_($i) [$ns_ node]	
	$node_($i) random-motion 0		;# disable random motion
}


#source $opt(nb)
# Define node movement model
puts "Loading connection pattern..."
source $opt(cp)
 
# Define traffic model
puts "Loading scenario file..."
source $opt(sc)

# Define node initial position in nam
#for {set i 0} {$i < $opt(nn)-3} {incr i} {
for {set i 0} {$i < $opt(nn)} {incr i} {
    # 20 defines the node size in nam, must adjust it according to your scenario
    # The function must be called after mobility model is defined
   $ns_ initial_node_pos $node_($i) 20
}

# Tell nodes when the simulation ends
for {set i 0} {$i < $opt(nn)} {incr i} {
    $ns_ at $opt(stop).000000001 "$node_($i) reset";
}

# tell nam the simulation stop time
$ns_ at  $opt(stop)	"$ns_ nam-end-wireless $opt(stop)"
$ns_ at  $opt(stop).000000001 "puts \"NS EXITING...\" ; $ns_ halt"
puts "Starting Simulation..."
$ns_ run
