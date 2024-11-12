#version 130

// TODO: add "in" variables
in vec3 eNormal;
in vec3 viewPos;

out vec4 fragColor;

// TODO: add uniform variables
uniform vec3 eLightDir;
uniform vec3 lightColor;
uniform float shininess; 
uniform vec3 diffuseCoeff;
uniform vec3 ambient;

void main()
{
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

    // Final color
    vec3 finalColor = ambient + diffuse + specular;
    fragColor = vec4(finalColor, 1.0);

}
