/*
**    TP CPE Lyon
**    Copyright (C) 2015 Damien Rohmer
**
**    This program is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#ifndef PLANE_HPP
#define PLANE_HPP

#include "primitive_basic.hpp"
#include "lib/3d/vec3.hpp"

namespace cpe
{

class plane : public primitive_basic
{
public:

    plane(vec3 const& position_param,vec3 const& normal_param);

    /** One point of the plane */
    vec3 const& position() const;
    /** Normal of the plance */
    vec3 const& normal() const;

    /** Intersection computation with a ray */
    bool intersect(ray const& ray_param,intersection_data& intersection) const override;


private:

    /** One point belonging to the plane */
    vec3 position_data;
    /** Normal of the plane */
    vec3 normal_data;
};

}

#endif
