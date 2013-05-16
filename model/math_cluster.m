coords = Import["/Users/hayssam/Documents/ISOP_0.2/tests/prior_coords.tsv"];
clusterables = (#[[3 ;;]] -> #[[;; 2]]) & /@ coords;
clusts = FindClusters[clusterables,DistanceFunction -> CosineDistance];
(*
(* Force at least two clusters *)
clusts =
  If[Length[clusts] == 1, 
   FindClusters[clusterables, 2, DistanceFunction -> CosineDistance],
   clusts];
  *)
Export["clusters.tsv", clusts[[All, All, 2]]]