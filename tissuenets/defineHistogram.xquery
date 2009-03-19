
declare variable $allnodes as xs:string external;
(: Farsight Input file that contains objects information :)

declare variable $objecttype as xs:string external;
(: Object type can be one of these:                       :) 
(: Microglia Endothelials Neurons Endothelials Astrocytes :)

declare variable $featurename as xs:string external;
(: Sample Feature: volume Texture                         :)


element histogram_data {
        attribute Class_Membership {$objecttype},
        attribute feature {$featurename},        

let $NCFeatures := doc($allnodes)/Farsight_Output/Nuclei/Nuclear_Features[@Class_Membership = $objecttype]
for $attr in $NCFeatures//@*
where name($attr)=$featurename
return <d>{data($attr)}</d>
}