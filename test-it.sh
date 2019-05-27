#!/bin/bash

cmd_warmup="julia warmup.jl"
eval $cmd_warmup&>warmup.log

export InFile1="test-dataset/scenario_1/case.con"
export InFile2="test-dataset/scenario_1/case.inl"
export InFile3="test-dataset/scenario_1/case.raw"
export InFile4="test-dataset/scenario_1/case.rop"
export NetworkModel="IEEE 14"

cp test-dataset/scenario_1/solution1.txt .

cmd="Julia -e 'include(\"MyJulia2.jl\"); MyJulia2(\"${InFile1}\", \"${InFile2}\", \"${InFile3}\", \"${InFile4}\", 34200, 2, \"${NetworkModel}\")'"
eval $cmd&>MyJulia2.log

rm solution1.txt
