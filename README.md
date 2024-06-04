# AOMDV-CLUSTER: An Energy-Efficient Clustering Algorithm for Wireless Sensor Networks

This repository contains the source code for `AOMDV-CLUSTER`, a novel clustering algorithm for wireless sensor networks (WSNs) that optimizes energy consumption and extends network lifetime.

## Algorithm Description

`AOMDV-CLUSTER` operates in two main steps:

1. **Configuration Step**: The base station creates clusters based on Euclidean distance and a threshold `K` (default value is 250, but it can be modified in the code). 

2. **Communication Step**: Once the clusters are created, the communication starts. The base station (BS) communicates with the cluster heads (CHs), CHs communicate with each other and the BS, and nodes communicate only with their own CH.

The algorithm is implemented by modifying the source code of the Ad hoc On-Demand Multipath Distance Vector (AOMDV) routing protocol in NS2.35. The energy model is defined and the configuration step is added in the `sendHELLO` function. The condition for starting the communication is added in the `recvRequest` function.

## How to Run

To run the simulation with this code, follow these steps:

1. Copy the `AOMDV-CLUSTER.cc` file and place it in the `/YOUR-NAME/ns-allinone-2.35/ns-2.35/aomdv.cc` directory.
2. Define the vector data in line 2029 and declare the number of nodes in line 121.
3. Compile the code.
4. You will need the `aomdv.tcl` file, the communication file, and the scene file to run the simulation.
5. Launch the simulation.

The program will print the energy consumed at 90 seconds and at the end of the simulation. Whenever a CH dies, it will be reconfigured and updated to a member of the CH whose energy is the maximum.

## Contact

If you have any questions or need further assistance, feel free to open an issue in this repository.
