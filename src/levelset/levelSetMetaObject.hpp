/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

# ifndef __BITPIT_LEVELSET_META_HPP__
# define __BITPIT_LEVELSET_META_HPP__

namespace bitpit{

class LevelSetObject;

class LevelSetMetaObject : public LevelSetObject{
    public:
    LevelSetMetaObject(int);
    virtual ~LevelSetMetaObject();

    bool            isPrimary() const override;
    virtual int     getPrimaryObjectId( long ) const =0;
    virtual std::vector<const LevelSetObject *> getPrimaryObjects() const =0;

};



}

#endif 
