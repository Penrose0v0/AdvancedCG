#pragma once

#include "Material.h"

class PerfectSpecularMaterial : public Material
{
private:
	PerfectSpecularMaterial()
		: ks(1,1,1)
	{
	}

	PerfectSpecularMaterial(const PerfectSpecularMaterial &m)
		: ks(m.ks), Material(m)
	{
	}

public:
	static PerfectSpecularMaterial *CreateMaterial()
	{
		PerfectSpecularMaterial *m = new PerfectSpecularMaterial();
		material_cache.push_back( m );
		return m;
	}

	static PerfectSpecularMaterial *CreateMaterial(const PerfectSpecularMaterial *_m)
	{
		PerfectSpecularMaterial *m = new PerfectSpecularMaterial(*_m);
		material_cache.push_back( m );
		return m;
	}

	~PerfectSpecularMaterial()
	{
	}

	Material_Type getMaterialType() const { return Perfect_Specular_Type; }

	void clone(Material **m) const { *m = new PerfectSpecularMaterial(*this); }

	vec3 getSpecularCoeff() const { return ks; }
	void setSpecularCoeff(const vec3 &k) { ks = k; }
	void setSpecularCoeff(float r, float g, float b) { ks = vec3(r,g,b); }

private:
	vec3 ks;	// specular

};
