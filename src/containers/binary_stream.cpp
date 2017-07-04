/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include <cassert>

#include "binary_stream.hpp"

/*!
* Stream a string from the binary stream.
*
* \param[in] stream input stream
* \param[in] value is the string to be streamed
*/
template<>
bitpit::IBinaryStream& operator>>(bitpit::IBinaryStream &stream, std::string &value)
{
    int size = 0;
    stream.read(size);

    if (size > 0) {
        std::vector<char> buffer(size);
        stream.read(buffer);
        value.assign(buffer.data(), size);
    }

    return stream;
}

/*!
* Stream a string to the binary stream.
*
* \param[in] stream is the output stream
* \param[in] value is the string to be streamed
*/
template<>
bitpit::OBinaryStream& operator<<(bitpit::OBinaryStream &stream, const std::string &value)
{
    int size = value.size();
    stream.write(size);

    if (size > 0) {
        stream.write(value.c_str(), size);
    }

    return stream;
}


namespace bitpit {

/*!
* \class BinaryStream
* \ingroup BinaryStream
*
* \brief Base class for defining input and output binary streams.
*/

/*!
* Constructor.
*
* Initialize an empty binary stream.
*/
BinaryStream::BinaryStream()
    : m_size(0), m_pos(0)
{
}

/*!
* Constructor.
*
* Initialize an empty binary stream with the assigned size.
*
* \param[in] size is the requested stream size
*/
BinaryStream::BinaryStream(std::size_t size)
    : m_pos(0)
{
    setSize(size);
}

/*!
* Constructor.
*
* Initialize a binary stream with the data contained in the specified buffer.
* The data is copied from the input buffer to the internal buffer.
*
* \param[in] buffer is the buffer that contains the data
* \param[in] capacity is the data capacity
*/
BinaryStream::BinaryStream(const char *buffer, std::size_t capacity)
{
    open(buffer, capacity);
}

/*!
* Constructor.
*
* Initialize a binary stream with the data contained in the specified buffer.
* The data is copied from the input buffer to the internal buffer.
*
* \param[in] buffer is the buffer that contains the data
*/
BinaryStream::BinaryStream(const std::vector<char> &buffer)
{
    open(buffer.data(), buffer.size());
}

/*!
* Open a binary stream initializing it with the data contained in the
* specified buffer.
*
* \param[in] buffer is the buffer that contains the data
* \param[in] size is the size (in bytes) of the data
*/
void BinaryStream::open(const char *buffer, std::size_t size)
{
    open(size);

    m_buffer.assign(buffer, buffer + size);
}

/*!
* Open a binary stream with the specified size.
*
* \param[in] size is the buffer size
*/
void BinaryStream::open(std::size_t size)
{
    m_pos = 0;

    setSize(size);
}

/*!
* Returns the cursor position within the stream.
*
* \result Returns the cursor position within the stream.
*/
std::streampos BinaryStream::tellg() const
{
    return m_pos;
}

/*!
* Set the cursor position within the stream.
*
* \param[in] pos is the cursor new position
* \result Returns true if new position is valid, false otherwise.
*/
bool BinaryStream::seekg(std::size_t pos)
{
    if (pos >= getSize()) {
        return false;
    }

    m_pos = pos;

    return true;
}

/*!
* Get a pointer the internal data.
*
* \return A pointer to the internal data.
*/
char * BinaryStream::data()
{
    return m_buffer.data();
}

/*!
* Get a cnstant pointer the internal data.
*
* \return A constant pointer to the internal data.
*/
const char * BinaryStream::data() const
{
    return m_buffer.data();
}

/*!
* Returns true if end of file condition is met.
*
* \result Returns true if end of file condition is met, false otherwise.
*/
bool BinaryStream::eof() const
{
    return (m_pos >= getSize());
}

/*!
* Set the cursor position within the stream.
*
* \param[in] offset is the offset value, relative to the way parameter
* \param[in] way is the offset direction
* \result Returns rue if new position is valid, false otherwise.
*/
bool BinaryStream::seekg(std::streamoff offset, std::ios_base::seekdir way)
{
    if ((way == std::ios_base::beg) && (offset < (long) getSize())) {
        m_pos = offset;
    } else if ((way == std::ios_base::cur) && (m_pos + offset < getSize())) {
        m_pos += offset;
    } else if ((way == std::ios_base::end) && ((long) getSize() - offset >= 0)) {
        m_pos = getSize() - offset;
    } else {
        return false;
    }

    return true;
}

/*!
* Return the size of the stream, expressed in bytes.
*
* \result The size of the stream, expressed in bytes.
*/
std::size_t BinaryStream::getSize() const
{
    return m_size;
}

/*!
* Set the size of the stream, expressed in bytes.
*
* \param[in] size is the new size (in bytes) of the stream
*/
void BinaryStream::setSize(std::size_t size)
{
    setCapacity(size);

    m_size = size;
}

/*!
* Requests that the stream capacity be at least the specified number of bytes.
*
* \param[in] capacity is the new size (in bytes) of the stream
*/
void BinaryStream::setCapacity(std::size_t capacity)
{
    m_buffer.resize(capacity);
}

/*!
* Returns the size of the storage space currently allocated for the stream,
* expressed in bytes.
*
* Capacity is not necessarily equal to the stream size; it can be equal
* or greater. The extra space can not be used for storing data, but
* can be used to guarantee a minimum internal storage size.
*
* Notice that this capacity does not suppose a limit on the size of the stream.
* When this capacity is exhausted and more is needed, it is automatically
* expanded by the (reallocating it storage space).
*
* \return The size of the storage space currently allocated for the stream,
* expressed in bytes.
*/
std::size_t BinaryStream::getCapacity() const
{
    return m_buffer.size();
}

/*!
* \class IBinaryStream
* \ingroup BinaryStream
*
* \brief Output binary stream.
*/

/*!
* Constructor.
*
* Initialize an empty binary stream.
*/
IBinaryStream::IBinaryStream()
    : BinaryStream()
{
}

/*!
* Constructor.
*
* Initialize an empty binary stream with the assigned size.
*
* \param[in] size is the stream size
*/
IBinaryStream::IBinaryStream(std::size_t size)
    : BinaryStream(size)
{
}

/*!
* Constructor.
*
* Initialize a binary stream with the data contained in the specified buffer.
* The data is copied from the input buffer to the internal buffer.
*
* \param[in] buffer is the buffer that contains the data
* \param[in] size is the data size
*/
IBinaryStream::IBinaryStream(const char *buffer, std::size_t size)
    : BinaryStream(buffer, size)
{
}

/*!
* Constructor.
*
* Initialize a binary stream with the data contained in the specified buffer.
* The data is copied from the input buffer to the internal buffer.
*
* \param[in] buffer is the buffer that contains the data
*/
IBinaryStream::IBinaryStream(const std::vector<char> &buffer)
    : BinaryStream(buffer)
{
}

/*!
* Open a binary stream with the specified size.
*
* \param[in] size is the buffer size
*/
void IBinaryStream::open(std::size_t size)
{
    BinaryStream::open(size);
}

/*!
* Open a binary stream initializing it with the data contained in the
* specified buffer.
*
* \param[in] buffer is the buffer that contains the data
* \param[in] size is the buffer size
*/
void IBinaryStream::open(const char *buffer, std::size_t size)
{
    BinaryStream::open(buffer, size);
}

/*!
* Read data from the stream.
*
* \param[out] data is the memory location that will contain the data
* \param[in] size is the size (in bytes) of the data to be read from the
* stream
*/
void IBinaryStream::read(char *data, std::size_t size)
{
    if (eof() || (m_pos + size) > getSize()) {
        throw std::runtime_error("Bad memory access!");
    }

    std::memcpy(reinterpret_cast<void*>(data), m_buffer.data() + m_pos, size);
    m_pos += size;
}

/*!
* \class OBinaryStream
* \ingroup BinaryStream
*
* \brief Output binary stream.
*/

/*!
* Constructor.
*
* Initialize an empty binary stream.
*/
OBinaryStream::OBinaryStream()
    : BinaryStream(), m_expandable(true)
{
}

/*!
* Constructor.
*
* Initialize a binary stream with the specified size.
*
* \param[in] capacity is the size of the stream
*/
OBinaryStream::OBinaryStream(std::size_t size)
    : BinaryStream(size), m_expandable(size == 0)
{
}

/*!
* Open a binary stream with the specified size.
*
* \param[in] capacity is the buffer size
*/
void OBinaryStream::open(std::size_t size)
{
    // Set the stream expandable to allow changing its capacity
    m_expandable = true;

    // Open a stream with the desidered capacity
    BinaryStream::open(size);

    // Set the definitive expandable flag
    m_expandable = (size == 0);
}

/*!
* Set the size of the stream, expressed in bytes.
*
* If the buffer is not expandable an exception is thrown.
*
* \param[in] size is the new size (in bytes) of the stream
*/
void OBinaryStream::setSize(std::size_t size)
{
    if (!m_expandable) {
        throw std::runtime_error("Stream is not expandable.");
    }

    BinaryStream::setSize(size);
}

/*!
* Requests the stream to reduce its size to fit the data currently contained
* in the stream.
*/
void OBinaryStream::squeeze()
{
    bool expandableFlag = m_expandable;

    // Set the stream expandable to allow changing its size
    m_expandable = true;

    // Update the size of the stream
    setSize(tellg());

    // Set the expandable flag back to its initial value
    m_expandable = expandableFlag;
}

/*!
* Write data into the stream.
*
* \param[in] data is the memory that contain the data
* \param[in] size is the size (in bytes) of the data to be written into
* the stream
*/
void OBinaryStream::write(const char *data, std::size_t size)
{
    if (getSize() - m_pos < size) {
        // If the stream is not expandable, the request for a new size
        // will throw an exception.
        std::size_t bufferSize = getSize() + size;
        setSize(bufferSize);
    }

    for (std::size_t i = 0;  i < size; ++i) {
        m_buffer[m_pos] = data[i];
        ++m_pos;
    }
}

}
