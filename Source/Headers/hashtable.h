//from https://github.com/preshing/CompareIntegerMaps for testing
#pragma once


#define FIRST_CELL(hash) (m_cells + ((hash) & (m_arraySize - 1)))
#define CIRCULAR_NEXT(c) ((c) + 1 != m_cells + m_arraySize ? (c) + 1 : m_cells)
#define CIRCULAR_OFFSET(a, b) ((b) >= (a) ? (b) - (a) : m_arraySize + (b) - (a))


inline uint32_t upper_power_of_two(uint32_t v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

inline uint64_t upper_power_of_two(uint64_t v)
{
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
}

// from code.google.com/p/smhasher/wiki/MurmurHash3
inline uint32_t integerHash(uint32_t h)
{
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

// from code.google.com/p/smhasher/wiki/MurmurHash3
inline uint64_t integerHash(uint64_t k)
{
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    return k;
}

//----------------------------------------------
//  HashTable
//
//  Maps pointer-sized integers to pointer-sized integers.
//  Uses open addressing with linear probing.
//  In the m_cells array, key = 0 is reserved to indicate an unused cell.
//  Actual value for key 0 (if any) is stored in m_zeroCell.
//  The hash table automatically doubles in size when it becomes 75% full.
//  The hash table never shrinks in size, even after Clear(), unless you explicitly call Compact().
//----------------------------------------------
class HashTable
{
public:
    typedef std::pair<unsigned, unsigned> Cell;
    /*
    struct Cell
    {
        unsigned key;
        unsigned value;
    };
    */
private:
    Cell* m_cells;
    size_t m_arraySize;
    size_t m_population;
    bool m_zeroUsed;
    Cell m_zeroCell;
    
    void Repopulate(size_t desiredSize)
    {
        assert((desiredSize & (desiredSize - 1)) == 0);   // Must be a power of 2
        assert(m_population * 4 <= desiredSize * 3);
        Cell* oldCells = m_cells;
        Cell* end = m_cells + m_arraySize;
        m_arraySize = desiredSize;
        m_cells = new Cell[m_arraySize];
        memset(m_cells, 0, sizeof(Cell) * m_arraySize);
        for (Cell* c = oldCells; c != end; c++)
        {
            if (c->first)
            {
                for (Cell* cell = FIRST_CELL(integerHash(c->first));; cell = CIRCULAR_NEXT(cell))
                {
                    if (!cell->first)
                    {
                        *cell = *c;
                        break;
                    }
                }
            }
        }
        delete[] oldCells;
    }


public:
    HashTable(size_t initialSize = 512)
    {
        // Initialize regular cells
        m_arraySize = initialSize;
        assert((m_arraySize & (m_arraySize - 1)) == 0);   // Must be a power of 2
        m_cells = new Cell[m_arraySize];
        memset(m_cells, 0, sizeof(Cell) * m_arraySize);
        m_population = 0;

        // Initialize zero cell
        m_zeroUsed = 0;
        m_zeroCell.first = 0;
        m_zeroCell.second = 0;
    }
    ~HashTable()
    {
        // Delete regular cells
        delete[] m_cells;
    }

    // Basic operations
    
    
    const Cell* Lookup(unsigned key) const
    {
        if (key)
        {
            // Check regular cells
            for (Cell* cell = FIRST_CELL(integerHash(key));; cell = CIRCULAR_NEXT(cell))
            {
                if (cell->first == key)
                    return cell;
                if (!cell->first)
                    return NULL;
            }
        }
        else
        {
            // Check zero cell
            if (m_zeroUsed)
                return &m_zeroCell;
            return NULL;
        }
    }
    
    Cell* Lookup(unsigned key)
    {
        if (key)
        {
            // Check regular cells
            for (Cell* cell = FIRST_CELL(integerHash(key));; cell = CIRCULAR_NEXT(cell))
            {
                if (cell->first == key)
                    return cell;
                if (!cell->first)
                    return NULL;
            }
        }
        else
        {
            // Check zero cell
            if (m_zeroUsed)
                return &m_zeroCell;
            return NULL;
        }
    };

    Cell* Insert(unsigned key)
    {
        if (key)
        {
            // Check regular cells
            for (;;)
            {
                for (Cell* cell = FIRST_CELL(integerHash(key));; cell = CIRCULAR_NEXT(cell))
                {
                    if (cell->first == key)
                        return cell;        // Found
                    if (cell->first == 0)
                    {
                        // Insert here
                        if ((m_population + 1) * 4 >= m_arraySize * 3)
                        {
                            // Time to resize
                            Repopulate(m_arraySize * 2);
                            break;
                        }
                        ++m_population;
                        cell->first = key;
                        return cell;
                    }
                }
            }
        }
        else
        {
            // Check zero cell
            if (!m_zeroUsed)
            {
                // Insert here
                m_zeroUsed = true;
                if (++m_population * 4 >= m_arraySize * 3)
                {
                    // Even though we didn't use a regular slot, let's keep the sizing rules consistent
                    Repopulate(m_arraySize * 2);
                }
            }
            return &m_zeroCell;
        }
    }

    void Delete(Cell* cell)
    {
        if (cell != &m_zeroCell)
        {
            // Delete from regular cells
            assert(cell >= m_cells && cell - m_cells < m_arraySize);
            assert(cell->first);

            // Remove this cell by shuffling neighboring cells so there are no gaps in anyone's probe chain
            for (Cell* neighbor = CIRCULAR_NEXT(cell);; neighbor = CIRCULAR_NEXT(neighbor))
            {
                if (!neighbor->first)
                {
                    // There's nobody to swap with. Go ahead and clear this cell, then return
                    cell->first = 0;
                    cell->second = 0;
                    m_population--;
                    return;
                }
                Cell* ideal = FIRST_CELL(integerHash(neighbor->first));
                if (CIRCULAR_OFFSET(ideal, cell) < CIRCULAR_OFFSET(ideal, neighbor))
                {
                    // Swap with neighbor, then make neighbor the new cell to remove.
                    *cell = *neighbor;
                    cell = neighbor;
                }
            }
        }
        else
        {
            // Delete zero cell
            assert(m_zeroUsed);
            m_zeroUsed = false;
            cell->second = 0;
            m_population--;
            return;
        }
    }

    void Clear()
    {
        // (Does not resize the array)
        // Clear regular cells
        memset(m_cells, 0, sizeof(Cell) * m_arraySize);
        m_population = 0;
        // Clear zero cell
        m_zeroUsed = false;
        m_zeroCell.second = 0;
    }


    void Compact()
    {
        Repopulate(upper_power_of_two((m_population * 4 + 3) / 3));
    }

    void Delete(unsigned key)
    {
        Cell* value = Lookup(key);
        if (value)
            Delete(value);
    }

    //----------------------------------------------
    //  Iterator
    //----------------------------------------------
    friend class Iterator;
    template<typename T>
    class Iterator
    {
    private:
        const HashTable* m_table;
        const T* m_cur;

    public:
        Iterator(const Iterator& b) : m_table(b.m_table), m_cur(b.m_cur) {}
        Iterator& operator =(Iterator const& b)
        {
            m_table = b.m_table;
            m_cur = b.m_cur;
            return *this;
        }
        Iterator(const HashTable &table) : m_table(&table)
        {
            m_cur = &m_table->m_zeroCell;
            if (!m_table->m_zeroUsed)
                m_cur = Next();
        }
        Iterator(const HashTable& table, const T* cur) : m_table(&table), m_cur(cur) {}
        inline Iterator& operator ++() { m_cur = Next(); return *this; }
        inline Iterator operator ++(int) { const T* old = m_cur; m_cur = Next(); return Iterator(*m_table, old); }
        const T* Next() const
        {
            const T* m_ret = m_cur;
            if (!m_ret) return m_ret;
            if (m_ret == &m_table->m_zeroCell)
                m_ret = &m_table->m_cells[-1];
            const T* end = m_table->m_cells + m_table->m_arraySize;
            while (++m_ret != end)
                if (m_ret->first) return m_ret;
            return m_ret = NULL;
        }

        inline const T* operator*() const { return m_cur; }
        inline const T* operator->() const { return m_cur; }
        inline bool operator ==(const Iterator& other) const { return m_cur == other.m_cur; }
        inline bool operator !=(const Iterator& other) const { return m_cur != other.m_cur; }
    };
    typedef typename Iterator<Cell> iterator;
    typedef typename Iterator<const Cell> const_iterator;

    inline iterator begin() { return iterator(*this); }
    inline iterator end() { return iterator(*this, NULL); }
    inline iterator find(unsigned key) { return iterator(*this, Lookup(key)); }
    inline const_iterator begin() const { return const_iterator(*this); }
    inline const_iterator end() const { return const_iterator(*this, NULL); }
    inline const_iterator find(unsigned key) const { return const_iterator(*this, Lookup(key)); }
};
