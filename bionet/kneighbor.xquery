
(:**************************************************************:)
(: This program translates an XML file to XGMML                 :)
(: Node definitions should be given with $xmlinput -XML file    :)
(: Edge definitions should be passed with $edges   -XML file    :)
(:**************************************************************:)

declare namespace math="java:java.lang.Math";

declare variable $nodes1 as xs:string external;
(: Farsight Input file that contains objects information :)

declare variable $nodes2 as xs:string external;
(: Another node definition file :)

declare variable $n as xs:integer external;
(: Number of n nearest neighbours :)

(:**************************************************************:)
(: computeDistances first finds distances from each object to   :) 
(: every other object given by the input file. Then these links :)
(: constructed from object-to-object distances are sorted       :)
(: The result is an XML element that contains all links and     :)
(: distances.                                                   :)
(: $xmlinput is the name of the file that contains node ids and :)
(: their coordinates. Its format:                               :)
(: <objects>                                                    :)
(:   <Nuclear_Features ID="3" X="570" Y="1020" Z="28" ...       :) 
(:                                                              :)
(:**************************************************************:)

declare function local:computeDistances($nodes1 as xs:string, $nodes2 as xs:string)
	as item()* {

element distances {
	attribute sourceclass {data(doc($nodes1)/objects/@class)},
	attribute targetclass {data(doc($nodes2)/objects/@class)},
for $currObject in doc($nodes1)/objects/Nuclear_Features
	return
     for $otherObject in doc($nodes2)/objects/Nuclear_Features
		let $dx := $otherObject/@X -  $currObject/@X 
		let $dy := $otherObject/@Y -  $currObject/@Y 
		let $dz := $otherObject/@Z -  $currObject/@Z 
                
                (: The voxel resolution for all images in FARSIGHT        :)
                (: is 0.36um X 0.36um X 1.5 um along the x, y, and z axes :)
                (: This affects distance computations                     :)

               	let $dxx := $dx * $dx * 0.36 * 0.36
               	let $dyy := $dy * $dy * 0.36 *0.36
               	let $dzz := $dz * $dz * 1.5 * 1.5

(: sqrt computing is expensive. We postbone doing it!      :)
(:             	let $distance := math:sqrt($dxx+$dyy+$dzz) :)
             	let $distance :=$dxx+$dyy+$dzz
	     where ($otherObject/@X != $currObject/@X)
	     order by $distance	
	     return        
                element edge {
		        element ID1 {data($currObject/@ID)},
		        element ID2 {data($otherObject/@ID)},
		        element X {data($currObject/@X)},
		        element Y {data($currObject/@Y)},
		        element Z {data($currObject/@Z)},
                        element dist {$distance}
		}
}
};

if (($nodes1 !='') and ($nodes2 != '')) then 
let $sortedDistances := local:computeDistances($nodes1, $nodes2)
return
element distances {
	attribute sources {data(doc($nodes1)/objects/@class)},
	attribute targets {data(doc($nodes2)/objects/@class)},

  for $currObject in 
	distinct-values($sortedDistances/descendant-or-self::distances/edge/ID1)
	let $i := $sortedDistances/descendant-or-self::distances/edge[ID1 = $currObject]
	order by ($currObject cast as xs:integer)
	return 
		for $j in (1 to $n)
		return
			element edge {
		        	element ID1 {data($i[$j]/ID1)},
		        	element ID2 {data($i[$j]/ID2)},
                                (: Convert Voxels to um's :)
		        	element X {round(0.36 * data($i[$j]/X))},
		        	element Y {round(0.36 * data($i[$j]/Y))},
		        	element Z {round(1.5 * data($i[$j]/Z))},
				(: We now can compute the real distance with sqrt :)
                        	element dist {math:sqrt(data($i[$j]/dist))}
		}
(:			subsequence($i,$j,1) :)

}
else
<distances> Node data is missing</distances>