#! /usr/bin/env bash

#root -l -b -q runCascAnalysis.C\(10,40,70,100,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(40,60,60,100,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(60,70,40,100,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(70,80,40,80,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(80,100,40,70,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(0,20,40,60,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(10,30,30,70,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(20,40,30,50,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(30,50,30,50,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(50,100,0,30,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(0,10,20,30,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(10,20,10,30,\"XiMinus\",\"SPDClusters\",\"V0M\"\)--> DCACascDaughters-3
#root -l -b -q runCascAnalysis.C\(20,30,0,20,\"XiMinus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(30,50,0,10,\"XiMinus\",\"SPDClusters\",\"V0M\"\) --> CascRadius-1

#root -l -b -q runCascAnalysis.C\(10,40,70,100,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(40,60,60,100,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(60,70,40,100,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(70,80,40,80,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(80,100,40,70,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(0,20,40,60,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(10,30,30,70,\"XiPlus\",\"SPDClusters\",\"V0M\"\) -->V0Mass 4
#root -l -b -q runCascAnalysis.C\(20,40,30,50,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(30,50,30,50,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
#root -l -b -q runCascAnalysis.C\(50,100,0,30,\"XiPlus\",\"SPDClusters\",\"V0M\"\) --> DCAV0ToPV-4
root -l -b -q runCascAnalysis.C\(0,10,20,30,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
root -l -b -q runCascAnalysis.C\(10,20,10,30,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
root -l -b -q runCascAnalysis.C\(20,30,0,20,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
root -l -b -q runCascAnalysis.C\(30,50,0,10,\"XiPlus\",\"SPDClusters\",\"V0M\"\)
