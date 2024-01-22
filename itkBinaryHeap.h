#ifndef itkBinaryHeap_h
#define itkBinaryHeap_h

#include "itkIndex.h"
#include "itkLevelSetNode.h"
#include "itkMacro.h"
#include "itkLightObject.h"
#include "itkIndex.h"
#include "itkNumericTraits.h"

#include <vector>

namespace itk
{
  /**
   * \class BinaryHeap
   * \brief Implements a binary heap tree to use with Fast Marching methods
   * 
   * BinaryHeap is a binary heap of LeveSetNode objects, where the root node 
   * is the one with the least value. 
   * 
   * Be careful if it needs to be debugged. Take into account that the node
   * container is an array, where the 0 position is there so that 
   */
  template< typename TPixel, unsigned int Dimension = 2 >
  class BinaryHeap : public LightObject
  {
  public:
    /** Standard class typedefs */
    typedef BinaryHeap 			Self;
    typedef LightObject			Superclass;
    typedef SmartPointer<Self>		Pointer;
    typedef SmartPointer<const Self>	ConstPointer;
      
    /** Heap-related typedefs */
    typedef TPixel 				PixelType;
    typedef LevelSetNode<PixelType,Dimension>	NodeType;
    typedef typename NodeType::IndexType	IndexType;

    typedef std::vector< NodeType>		HeapContainer;
    
    /** Method for creation through the object factory. */
    itkNewMacro( Self );

    /** Run-time type information (and related methods). */
    itkTypeMacro( BinaryHeap, LightObject );
    

    /** Check whether the heap is empty*/
    bool IsEmpty() { return (m_CurrentSize==0); }
    
    /** Get the number of nodes in the heap */
    unsigned long Size() { return m_CurrentSize; }
    
    /** Get the node at the top */
    NodeType Top()
    {
      if (m_CurrentSize>0)
	return m_Heap[1];
      else
	itkExceptionMacro( "Empty Heap" );
    }
    
    /** Remove the node at the top. It does not return this node */
    void Pop()
    {
      //Put the last node on the top of the heap and bubble it down.
      m_Heap[1] = m_Heap[ m_CurrentSize ];
      --m_CurrentSize;
      
      this->BubbleDown( 1 );
      
    }
    
    /** Inserts node into the heap */
    void Insert( NodeType & x )
    {
      //Make sure there is enough capacity in the vector
      if (m_Heap.capacity() <= m_CurrentSize+1)
      {
	m_Heap.reserve( 2*m_CurrentSize );
      }
      
      //Add the node to the end of the heap vector and bubble it up
      m_Heap.insert(m_Heap.begin()+1+m_CurrentSize, x);
      ++m_CurrentSize;
      
      BubbleUp( m_CurrentSize );
      
    }
    
    /** If a node changes its value and it is smaller, bubble it up */
    void RefreshIfSmaller( NodeType & x )
    {
      //Find the node
      IndexType idx = x.GetIndex();
      
      unsigned long pos = 0;
      for (unsigned long i=1;i<=m_CurrentSize;i++)
      {
	if ( idx == m_Heap[i].GetIndex() )
	{
	  pos = i;
	  break;
	}
      }
      
      //Bubble it up if needed
      if ((pos>0)&&(x<m_Heap[pos]))
	BubbleUp(pos);
      else
	itkExceptionMacro( "Node not found." );
    }
    
    void ResetHeap()
    {
      //If there is memory allocated, it's better not to free it.
      m_Heap.resize(1); //Preserve the empty node at the beginning
      m_CurrentSize = 0;
    }
    
  protected:
    BinaryHeap()
    {
      NodeType vacio; 
      vacio.SetValue( NumericTraits<PixelType>::min() );
      
      m_Heap.clear();
      m_Heap.reserve( 100 );
      
      m_Heap.push_back( vacio );
      m_CurrentSize = 0;
    }
    virtual ~BinaryHeap() {}
    
    void BubbleDown( unsigned long hole )
    {
      unsigned long child;
      NodeType tmp = m_Heap[hole];
      
      for ( ; 2*hole<=m_CurrentSize ; hole=child )
      {
	child = 2*hole;
	//Check we are looking at the smallest child
	if ((child!=m_CurrentSize) && (m_Heap[child] > m_Heap[child+1]))
	  ++child;
	
	//If the child is smaller than the parent, swap them. 
	//If not, we cannot go any further.
	if (m_Heap[child] < tmp)
	  m_Heap[ hole ] = m_Heap[ child ];
	else
	  break;
      }
      m_Heap[hole] = tmp;
      
    }
    
    void BubbleUp( unsigned long hole )
    {
      unsigned long parent = hole/2;
      NodeType x = m_Heap[hole];
      
      while ( (parent>0) && (x < m_Heap[parent]))
      {
	m_Heap[hole] = m_Heap[parent];
	hole = parent;
	parent = hole/2;
      }
      m_Heap[hole] = x;
    }
    
  private:
    //Heap structure
    HeapContainer m_Heap;
    unsigned long m_CurrentSize;
    
  };
  
}

#endif
