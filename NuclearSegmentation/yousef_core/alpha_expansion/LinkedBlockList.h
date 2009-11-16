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

/* Singly Linked List of Blocks */
// This data structure should be used only for the GCoptimization class implementation
// because it lucks some important general functions for general list, like remove_item()
// The head block may be not full 
// For regular 2D grids, it's better to set GCLL_BLOCK_SIZE to 2
// For other graphs, it should be set to the average expected number of neighbors
// Data in linked list for the neighborhood system is allocated in blocks of size GCLL_BLOCK_SIZE 

#ifndef __LINKEDBLOCKLIST_H__
#define __LINKEDBLOCKLIST_H__

#define GCLL_BLOCK_SIZE 4  
// GCLL_BLOCKSIZE should "fit" into the type BlockType. That is 
// if GCLL_BLOCKSIZE is larger than 255 but smaller than largest short integer
// then  BlockType should be set to short
typedef char BlockType;

//The type of data stored in the linked list
typedef void * ListType;

class LinkedBlockList{

public: 
	void addFront(ListType item);
	inline bool isEmpty(){if (m_head == 0) return(true); else return(false);};
	inline LinkedBlockList(){m_head = 0; m_head_block_size = GCLL_BLOCK_SIZE;}; 
	~LinkedBlockList();

	// Next three functins are for the linked list traversal
	inline void setCursorFront(){m_cursor = m_head; m_cursor_ind = 0;};
	ListType next();
	bool hasNext();

	//added by yousef al-kofahi (RPI) on may 24th 2008
	inline int getBlockSize(){return (int)m_head_block_size;};

private:
	typedef struct LLBlockStruct{
		ListType m_item[GCLL_BLOCK_SIZE];
		struct LLBlockStruct *m_next;
	} LLBlock;

	LLBlock *m_head;
	// Remembers the number of elements in the head block, since it may not be full
	BlockType m_head_block_size;
	// For block traversal, points to current element in the current block
	BlockType m_cursor_ind;
	// For block traversal, points to current block in the linked list
	LLBlock *m_cursor;
};

#endif

