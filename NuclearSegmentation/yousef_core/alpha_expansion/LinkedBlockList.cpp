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

#include "LinkedBlockList.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include "mex.h"

/*********************************************************************/

void LinkedBlockList::addFront(ListType item) {

	if ( m_head_block_size == GCLL_BLOCK_SIZE )
	{
		LLBlock *tmp      = (LLBlock *) new LLBlock;
		if ( !tmp )
			std::cerr<<"GraphCut:LinkedBlockList:addFront->out of mem"<<std::endl;
            //mexErrMsgIdAndTxt("GraphCut:LinkedBlockList:addFront","out of mem"); // BAGON REPLACED{printf("\nOut of memory");exit(1);}
		tmp -> m_next     = m_head;
		m_head            = tmp;
		m_head_block_size = 0;
	}
	
	m_head ->m_item[(int)m_head_block_size] = item;
	m_head_block_size++;
}

/*********************************************************************/

ListType LinkedBlockList::next()
{
	ListType toReturn = m_cursor -> m_item[(int)m_cursor_ind];

	m_cursor_ind++;

	if ( m_cursor == m_head && m_cursor_ind >= m_head_block_size )
	{
		m_cursor     = m_cursor ->m_next;
		m_cursor_ind = 0;
	}
	else if ( m_cursor_ind == GCLL_BLOCK_SIZE )
	{
		m_cursor = m_cursor ->m_next;
		m_cursor_ind = 0;
	}
	return(toReturn);
}

/*********************************************************************/

bool LinkedBlockList::hasNext()
{
	if ( m_cursor != 0 ) return (true);
	else return(false);
}


/*********************************************************************/

LinkedBlockList::~LinkedBlockList()
{
	LLBlock *tmp;

	while ( m_head != 0 ) 
	{
		tmp = m_head;
		m_head = m_head->m_next;
		delete tmp;
	}
};

/*********************************************************************/

