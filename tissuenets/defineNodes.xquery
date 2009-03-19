
declare variable $allnodes as xs:string external;
(: Farsight Input file that contains objects information :)

declare variable $objecttype as xs:string external;
(: Object type can be one of these:                       :) 
(: Microglia Endothelials Neurons Endothelials Astrocytes :)

<objects> {
	attribute class {$objecttype},
for $fobject in 
   doc($allnodes)/Farsight_Output/Nuclei/Nuclear_Features[@Class_Membership = $objecttype]
return  $fobject

}</objects>

