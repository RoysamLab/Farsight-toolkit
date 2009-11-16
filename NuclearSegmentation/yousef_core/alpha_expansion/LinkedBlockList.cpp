/* 
 * Copyright 2009 Rensselaer Polytechnic Institute
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but 
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 * for more details.
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the Free Software Foundation, Inc., 
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

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
		tmp = NULL;
	}
};

/*********************************************************************/

