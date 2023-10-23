#numnodes=80
#graphid=5015
numnodes=51
graphid=4577
drawgraph=1
subgraph=1
python3 stratisfimal.py $numnodes 1 $drawgraph 1 1 60 1 scripts/random_graph.txt 0 $graphid $subgraph
python3 stratisfimal.py $numnodes 1 $drawgraph 1 1 60 1 scripts/random_graph.txt 0 $graphid 0
