#version 130

in vec3 vVertexNormal;
out vec4 fragColor;

void main()
{
	// TODO: write an appropriate code here
	vec3 tmp = 0.5 * vVertexNormal + 0.5; 
	fragColor = vec4(tmp, 1.0); 

}
