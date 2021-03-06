/*
 * $Id$
 *
 * ***** BEGIN GPL LICENSE BLOCK *****
 *
 * Copyright 2009-2011 Jörg Hermann Müller
 *
 * This file is part of AudaSpace.
 *
 * Audaspace is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * AudaSpace is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Audaspace; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * ***** END GPL LICENSE BLOCK *****
 */

/** \file audaspace/intern/AUD_StreamBufferFactory.cpp
 *  \ingroup audaspaceintern
 */


#include "AUD_StreamBufferFactory.h"
#include "AUD_BufferReader.h"
#include "AUD_Buffer.h"

#include <cstring>

AUD_StreamBufferFactory::AUD_StreamBufferFactory(AUD_IFactory* factory) :
	m_buffer(new AUD_Buffer())
{
	AUD_IReader* reader = factory->createReader();

	m_specs = reader->getSpecs();

	int sample_size = AUD_SAMPLE_SIZE(m_specs);
	int length;
	int index = 0;
	sample_t* buffer;

	// get an approximated size if possible
	int size = reader->getLength();

	if(size <= 0)
		size = AUD_BUFFER_RESIZE_BYTES / sample_size;
	else
		size += m_specs.rate;

	// as long as we fill our buffer to the end
	while(index == m_buffer.get()->getSize() / sample_size)
	{
		// increase
		m_buffer.get()->resize(size*sample_size, true);

		// read more
		length = size-index;
		reader->read(length, buffer);
		memcpy(m_buffer.get()->getBuffer() + index * m_specs.channels,
			   buffer,
			   length * sample_size);
		size += AUD_BUFFER_RESIZE_BYTES / sample_size;
		index += length;
	}

	m_buffer.get()->resize(index * sample_size, true);
	delete reader;
}

AUD_IReader* AUD_StreamBufferFactory::createReader() const
{
	return new AUD_BufferReader(m_buffer, m_specs);
}
