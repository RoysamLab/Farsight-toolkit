
(:**************************************************************:)
(: This program translates an XML file to XGMML                 :)
(: Node definitions should be given with $xmlinput -XML file    :)
(: Edge definitions should be passed with $edges   -XML file    :)
(:**************************************************************:)

declare variable $nodes1 as xs:string external;
(: Farsight Input file that contains objects information :)

declare variable $nodes2 as xs:string external;
(: Farsight Input file that contains objects information :)

declare variable $links as xs:string external;
(: Edge definitions file :)

declare function local:createNodes($nodes as xs:string, $nodeshape as xs:string)
	as item()* {

	let $firstletter :=substring(data(doc($nodes)/objects/@class),1,1)
	for $fobject in doc($nodes)/objects/Nuclear_Features
	return  
	  element node {
		(:id should not change in farsight but for a better display :)
		(:the first letter of class type will be prefixed for now   :)
		(:attribute id {concat($firstletter,$fobject/@ID)},         :)

		attribute id {$fobject/@ID},

		(: label and weight should be decided :)
	  	attribute label {concat($firstletter,$fobject/@ID)},
	  	attribute weight {0},

	 	 element graphics {
			attribute type {$nodeshape},
                         (: Convert Voxels to um's :)
			attribute x {round(0.36 * $fobject/@X)},
			attribute y {round(0.36 * $fobject/@Y)},
			attribute z {round(1.5 * $fobject/@Z)}	
		}	
	}
};

element graph {
attribute directed {0},
attribute graphic {1},
attribute Layout {'points'},

(:**************************************************************:)
(: First create nodes                                           :)
(:**************************************************************:)

local:createNodes($nodes1,'circle'),

(: If there are two types of nodes, their shapes must be different :)
if ($nodes1 != $nodes2) then 
	local:createNodes($nodes2,'ver_ellipsis')
(: rhombus :)
else (),

(:**************************************************************:)
(: Create Links                                                 :)
(:**************************************************************:)

for $obj_to_obj in doc($links)/distances/edge
		return 
			element edge {
				attribute source {$obj_to_obj/ID1/text()},
	  			attribute target {$obj_to_obj/ID2/text()},
				attribute weight {$obj_to_obj/dist/text()}
			}
}