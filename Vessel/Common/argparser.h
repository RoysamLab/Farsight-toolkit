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

#ifndef __ARG_PARSER_H__
#define __ARG_PARSER_H__

#include "string.h"
#include "stdio.h"
#include "assert.h"
#include "stdlib.h"

enum RENDER_MODE { RENDER_MATERIALS, RENDER_LIGHTS, RENDER_SELECTION,RENDER_UNDISTRIBUTED, RENDER_ABSORBED, RENDER_RADIANCE, RENDER_FORM_FACTORS };
//#define NUM_RENDER_MODES 6
//enum RENDER_MODE { RENDER_MATERIALS, RENDER_LIGHTS, RENDER_UNDISTRIBUTED, RENDER_ABSORBED, RENDER_RADIANCE, RENDER_FORM_FACTORS };
#define NUM_RENDER_MODES 2


class ArgParser {

public:

  ArgParser() { DefaultValues(); }

  ArgParser(int argc, char *argv[]) {
    DefaultValues();

    for (int i = 1; i < argc; i++) {
      if (!strcmp(argv[i],"-input")) {
	i++; assert (i < argc); 
	input_file = argv[i];
      } else if (!strcmp(argv[i],"-size")) {
	i++; assert (i < argc); 
	width = height = atoi(argv[i]);
      } else if (!strcmp(argv[i],"-wireframe")) {
        wireframe = true;
      } else if (!strcmp(argv[i],"-num_samples")) {
	i++; assert (i < argc); 
	num_samples = atoi(argv[i]);
      } else if (!strcmp(argv[i],"-interpolate")) {
        interpolate = true;
      } else {
	printf ("whoops error with command line argument %d: '%s'\n",i,argv[i]);
	assert(0);
      }
    }
  }

  void DefaultValues() {
    input_file = NULL;
    width = 800;
    height = 800;
    wireframe = false;
    interpolate = false;
    render_mode = RENDER_RADIANCE;
    animate = false;
    num_samples = 1;
    tone_map = false;
    points = false;
    volume_rendering = false;
	edit_mode = false;
	show_number = false;
	votes = false;
	vote_number = 5;
  }

  // ==============
  // REPRESENTATION
  // all public! (no accessors)

  char *input_file;
  int width;
  int height;
  bool wireframe;
  bool interpolate;
  enum RENDER_MODE render_mode;
  bool animate;
  int num_samples;
  bool tone_map;
  bool points;
  bool volume_rendering;
  bool edit_mode;
  bool show_number;
  bool votes;
  int vote_number;
};

#endif
