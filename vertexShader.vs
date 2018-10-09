#version 330 core
out vec4 vertexColor;
layout (location = 0) in vec3 aPos; 
layout (location = 1) in float densitySource; /* source of density */
uniform mat4 transform;
uniform mat4 projection;

void main()
{
	gl_Position = projection*transform * vec4(aPos, 1.0);
    if (densitySource > 0)
    {
        vertexColor = vec4(densitySource, densitySource, densitySource, densitySource);
    }
    else
    {
        vertexColor = vec4(1, 0.0, 0.0, 1.0);
    }
}

