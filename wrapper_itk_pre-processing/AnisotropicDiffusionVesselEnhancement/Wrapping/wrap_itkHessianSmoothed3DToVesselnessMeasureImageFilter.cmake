WRAP_CLASS("itk::SymmetricSecondRankTensor")
	WRAP_TEMPLATE("${ITKM_D}3" "${ITKT_D}, 3" )
END_WRAP_CLASS()

WRAP_CLASS("itk::Image" POINTER)
	WRAP_TEMPLATE("SSRT${ITKM_D}33" "itk::SymmetricSecondRankTensor< ${ITKT_D}, 3 >, 3")
END_WRAP_CLASS()

WRAP_CLASS("itk::ImageSource" POINTER)	
	WRAP_TEMPLATE("ISSRT${ITKM_D}33" "itk::Image< itk::SymmetricSecondRankTensor< ${ITKT_D}, 3>, 3 >")
END_WRAP_CLASS()

WRAP_CLASS("itk::ImageToImageFilter" POINTER)
	FOREACH(d "3")
		FOREACH(t ${WRAP_ITK_REAL})
      		WRAP_TEMPLATE("${ITKM_I${t}${d}}ISSRT${ITKM_D}${d}${d}" "${ITKT_I${t}${d}}, itk::Image< itk::SymmetricSecondRankTensor< ${ITKT_D}, ${d} >, ${d} >")
			WRAP_TEMPLATE("ISSRT${ITKM_D}${d}${d}${ITKM_I${t}${d}}" "itk::Image< itk::SymmetricSecondRankTensor< ${ITKT_D}, ${d} >, ${d} >, ${ITKT_I${t}${d}}")
		ENDFOREACH()
	ENDFOREACH()
END_WRAP_CLASS()


WRAP_CLASS("itk::HessianSmoothed3DToVesselnessMeasureImageFilter" POINTER)
  FOREACH(t ${WRAP_ITK_REAL})
  WRAP_TEMPLATE("${ITKM_${t}}" "${ITKT_${t}}")
  ENDFOREACH(t)
END_WRAP_CLASS()


