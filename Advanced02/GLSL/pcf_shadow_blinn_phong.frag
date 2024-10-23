#version 130

// TODO: define "in" variables
in vec3 viewPos; 
in vec3 eNormal;
in vec4 shadowTexCor; 

out vec4 fragColor;

uniform sampler2DShadow shadowTex;
uniform vec2 texMapScale;

// TODO: define uniform variables
uniform vec3 eLightDir;
uniform vec3 lightColor;
uniform float shininess; 
uniform vec3 diffuseCoeff;
uniform vec3 ambient;

float offsetLookup(sampler2DShadow map, vec4 loc, vec2 offset)
{
	return textureProj(map, vec4(loc.xy + offset * texMapScale * loc.w, loc.z, loc.w));
}

float samplingLimit = 3.5;

void main()
{
	// HINT: The visibility (i.e., shadowed or not) can be fetched using the "offsetLookup" function
	//       in the double loops. The y and x loops ranges from -samplingLimit to samplingLimit
	//       with a stepsize 1.0. Finally, the summed visibility is divided by 64 to normalize.

	// TODO: rewirte this function

	vec3 norm = normalize(eNormal);
    vec3 lightDir = normalize(eLightDir);

    // Diffuse
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * diffuseCoeff * lightColor;

    // Specular
    vec3 viewDir = normalize(-viewPos);
    vec3 halfVec = normalize(eLightDir + viewDir); 
    float spec = pow(max(dot(halfVec, norm), 0.0), shininess);
    vec3 specular = spec * lightColor;

    float shadowFactor = 0.0;
    int sampleCount = 0;

    // Loop through the neighboring texels within the sampling range
    for (float y = -samplingLimit; y <= samplingLimit; y += 1.0)
    {
        for (float x = -samplingLimit; x <= samplingLimit; x += 1.0)
        {
            shadowFactor += offsetLookup(shadowTex, shadowTexCor, vec2(x, y));
            sampleCount++;
        }
    }

    shadowFactor /= float(sampleCount);

    // Final color
    vec3 finalColor = ambient + shadowFactor * (diffuse + specular);
    fragColor = vec4(finalColor, 1.0);
}
