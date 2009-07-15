/*=========================================================================
Copyright 2009 Rensselaer Polytechnic Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. 
=========================================================================*/

#include "farsight_object.h"

template<unsigned dim >
typename farsight_object<dim>::ImageConstPtrType
farsight_object<dim>::
mouse_click(LocationType const& points)
{
  assert(!"mouse_click not defined!");
  return NULL;
}

template< unsigned dim >
void
farsight_object<dim>::
set_image( std::string const& image_path,
           std::string const& image_name)
{
  image_name_ = image_name;
  image_path_ = image_path;
}

template< unsigned dim >
std::string const&
farsight_object<dim>::
get_image_name()
{
  return image_name_;
}

template< unsigned dim >
std::string const&
farsight_object<dim>::
get_image_path()
{
  return image_path_;
}

template< unsigned dim >
typename farsight_object<dim>::ImageConstPtrType
farsight_object<dim>::
add_object(LocationListType const& points)
{
  assert(!"add_object not defined!");
  return NULL;
}

template< unsigned dim >
typename farsight_object<dim>::ImageConstPtrType
farsight_object<dim>::
add_object(LocationType const& point)
{
  assert(!"add_object not defined!");
  return NULL;
}
 
template< unsigned dim >
typename farsight_object<dim>::ImageConstPtrType
farsight_object<dim>::
delete_object(LocationType const& point )
{
  assert(!"delete_object not defined!");
  return NULL;
}

template< unsigned dim >
typename farsight_object<dim>::ImageConstPtrType
farsight_object<dim>::
delete_object(LocationListType const& point )
{
  assert(!"delete_object not defined!");
  return NULL;
}

template< unsigned dim >
typename farsight_object<dim>::ImageConstPtrType
farsight_object<dim>::
merge_objects(LocationListType const& points)
{
  assert(!"merge_object not defined!");
  return NULL;
}

template< unsigned dim >
typename farsight_object<dim>::ImageConstPtrType
farsight_object<dim>::
split_object(LocationType const& points )
{
  assert(!"merge_object not defined!");
  return NULL;
}

template< unsigned dim >
typename farsight_object<dim>::ImageConstPtrType
farsight_object<dim>::
split_object(LocationListType const& points )
{
  assert(!"merge_object not defined!");
  return NULL;
}

template<unsigned dim >
void 
farsight_object<dim>::
set_parameters( std::vector< std::string > const& parameters)
{
  parameters_ = parameters;
}

template<unsigned dim >
void 
farsight_object<dim>::
set_parameters( std::string const& xml_filename)
{
  assert(!"set_parameters not defined to read from an xml file!");
}

template< unsigned dim >
void 
farsight_object<dim>::
set_default_parameters()
{
  assert(!"set_default_parameters not defined!");
}

template< unsigned dim > 
typename farsight_object<dim>::ImageConstPtrType 
farsight_object<dim>::
display()
{
  return label_image_.GetPointer();
}

/*
template< unsigned dim >
bool 
farsight_object<dim>::run()
{
  assert(!"run() not defined!");
  return false;
}

template< unsigned dim >
void 
farsight_object<dim>::
write_xml(std::string const& filename)
{
  assert(!"write_xml() not defined!");
}

template< unsigned dim >
void 
farsight_object<dim>::
read_xml(std::string const& filename)
{
  assert(!"read_xml() not defined!");
}
*/
