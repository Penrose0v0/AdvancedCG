#pragma once

#include "Material.h"
#include <algorithm>

class BlinnPhongMaterial : public Material
{
private:
	BlinnPhongMaterial()
		: shininess(32), ks(1,1,1), kd(1,1,1)
	{
	}

	BlinnPhongMaterial(const BlinnPhongMaterial &m)
		: shininess(m.shininess), ks(m.ks), kd(m.kd), Material(m)
	{
	}

public:
	static BlinnPhongMaterial *CreateMaterial()
	{
		BlinnPhongMaterial *m = new BlinnPhongMaterial();
		material_cache.push_back( m );
		return m;
	}

	static BlinnPhongMaterial *CreateMaterial(const BlinnPhongMaterial *_m)
	{
		BlinnPhongMaterial *m = new BlinnPhongMaterial(*_m);
		material_cache.push_back( m );
		return m;
	}

	~BlinnPhongMaterial()
	{
	}

	Material_Type getMaterialType() const { return Blinn_Phong_Type; }

	Real getShininess() const { return shininess; }
	void setShininess(Real s) { shininess = s; }

	vec3 getSpecularCoeff() const { return ks; }
	void setSpecularCoeff(const vec3 &k) { ks = k; }
	void setSpecularCoeff(Real r, Real g, Real b) { ks = vec3(r,g,b); }

	vec3 getDiffuseCoeff() const { return kd; }
	void setDiffuseCoeff(const vec3 &k) { kd = k; }
	void setDiffuseCoeff(Real r, Real g, Real b) { kd = vec3(r,g,b); }

private:
	Real shininess; // power of cosine lobe
	vec3 ks;	// phong
	vec3 kd;	// diffuse

};
