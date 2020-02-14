#!/bin/bash

#cmd_warmup="julia warmup.jl"
#eval $cmd_warmup&>warmup.log

export InFile1="../../test/data/scenario_01/case.con"
export InFile2="../../test/data/scenario_01/case.inl"
export InFile3="../../test/data/scenario_01/case.raw"
export InFile4="../../test/data/scenario_01/case.rop"
export NetworkModel="IEEE 14"

cmd_one="julia -e 'include(\"MyJulia1.jl\"); MyJulia1(\"${InFile1}\", \"${InFile2}\", \"${InFile3}\", \"${InFile4}\", 600, 2, \"${NetworkModel}\")'"
echo $cmd_one
eval $cmd_one&>MyJulia1.log

cmd_two="julia -e 'include(\"MyJulia2.jl\"); MyJulia2(\"${InFile1}\", \"${InFile2}\", \"${InFile3}\", \"${InFile4}\", 600, 2, \"${NetworkModel}\")'"
echo $cmd_two
eval $cmd_two&>MyJulia2.log
