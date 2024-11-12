#pragma once

#include "Material.h"
#include <algorithm>

class DiffuseMaterial : public Material
{
private:
	DiffuseMaterial()
		: kd(1,1,1)
	{
	}

	DiffuseMaterial(const DiffuseMaterial &m)
		: kd(m.kd), Material(m)
	{
	}

public:
	static DiffuseMaterial *CreateMaterial()
	{
		DiffuseMaterial *m = new DiffuseMaterial();
		material_cache.push_back( m );
		return m;
	}

	static DiffuseMaterial *CreateMaterial(const DiffuseMaterial *_m)
	{
		DiffuseMaterial *m = new DiffuseMaterial(*_m);
		material_cache.push_back( m );
		return m;
	}

	~DiffuseMaterial()
	{
	}

	Material_Type getMaterialType() const { return Diffuse_Type; }

	void clone(Material **m) const { *m = new DiffuseMaterial(*this); }

	vec3 getDiffuseCoeff() const { return kd; }
	void setDiffuseCoeff(const vec3 &k) { kd = k; }
	void setDiffuseCoeff(Real r, Real g, Real b) { kd = vec3(r,g,b); }

private:
	vec3 kd;

};
