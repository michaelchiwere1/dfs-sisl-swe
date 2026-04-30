#!/bin/bash

#PJM -g M101
#PJM -L rscgrp=cxm
###PJM -L rscgrp=dbg
#PJM -L vnode=1
#PJM -L vnode-core=40
#PJM -L vnode-mem=160Gi
###PJM -L vnode-mem=320Gi
#PJM -L elapse=00:10:00
#PJM -s
###PJM -j
###PJM -P exec-policy=simplex

export OMP_NUM_THREADS=40

#export OMP_STACKSIZE=1024m
export OMP_STACKSIZE=32000m

./sw_eul_sh
