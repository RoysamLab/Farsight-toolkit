SET( fs_sources
  farsight_object.h             farsight_object.txx
  fsc_channel_accessor.h        fsc_channel_accessor.txx
)
AUX_SOURCE_DIRECTORY( Templates fs_sources )


ADD_LIBRARY( fscommon ${fs_sources} )
TARGET_LINK_LIBRARIES( fscommon  ${ITK_LIBRARIES} )
