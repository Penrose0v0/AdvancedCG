#pragma once

#include "Material.h"

class SpecularRefractionMaterial : public Material
{
private:
	SpecularRefractionMaterial()
		: eta(1.33f), ks(1,1,1)
	{
	}

	SpecularRefractionMaterial(const SpecularRefractionMaterial &m)
		: eta(m.eta), ks(m.ks), Material(m)
	{
	}

public:
	static SpecularRefractionMaterial *CreateMaterial()
	{
		SpecularRefractionMaterial *m = new SpecularRefractionMaterial();
		material_cache.push_back( m );
		return m;
	}

	static SpecularRefractionMaterial *CreateMaterial(const SpecularRefractionMaterial *_m)
	{
		SpecularRefractionMaterial *m = new SpecularRefractionMaterial(*_m);
		material_cache.push_back( m );
		return m;
	}

	~SpecularRefractionMaterial()
	{
	}

	Material_Type getMaterialType() const { return Specular_Refraction_Type; }

	vec3 getSpecularCoeff() const { return ks; }
	void setSpecularCoeff(const vec3 &k) { ks = k; }
	void setSpecularCoeff(float r, float g, float b) { ks = vec3(r,g,b); }

	float getRefractionIndex() const { return eta; }
	void setRefractionIndex(float _eta) { eta = _eta; }

private:
	float eta;
	vec3 ks;	// refraction

};
