/*
Date : 22/08/18
Command to compile : g++ -std=c++11 main.cc Shader.cc glad.c -lglfw -lGL -lX11 -lpthread -lXrandr -lXi -ldl -o hello
g++ -std=c++11 src/main.cc src/Shader.cc src/glad.c src/solverUtils.cc -lglfw3 -lGLU -lGL -lX11 -lXxf86vm -lpthread -lXrandr -lXi -ldl -lXinerama -lXcursor -lrt -lm -o hello
*/

#include <glad/glad.h>	// glad manages hardware specific function pointers for openGL, gives address of extensions etc
#include <GLFW/glfw3.h>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>	

#include "../includes/Shader.h"	
#include "../includes/solverUtils.h"

using namespace std;

//macros 

// #define WINDOW_SIZE 600
//array value accessor
#define IX(i,j) ((i)+(N+2)*(j))
// //resolution
#define N 100
// //overall grid size
// #define size = (N+2)*(N+2)	

// static u[size], v[size], u_prev[size], v_p rev[size];
// static dens[size], dens_prev[size]; 

// grid of density sources
float *dens_prev;
float *dens;
float *u, *v, *u_prev, *v_prev;

int WINDOW_SIZE_WIDTH = 600;
int WINDOW_SIZE_HEIGHT = 600;

// vertices for a grid
float gridPoints[4*(N+2)*(N+2)];

#define size (N+2)*(N+2)
float vel_pos[6*(N+2)*(N+2)];

float dens_grid[4*size];

static float dt = 0.00004f;
static float source = 200.0f;
static float diff = 0.0004f;//0.0000000000000000000000000001f;

// callback function when window is resized to adjust viewport
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
	WINDOW_SIZE_WIDTH = width;
	WINDOW_SIZE_HEIGHT = height;
}  

void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
	{
        glfwSetWindowShouldClose(window, true);
	}
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
	{
		std::cout << "key e pressed" << std::endl;
	}
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	// std::cout << xpos << ", " << ypos  << std::endl;
}

// remap values from mouse coordinates to screen/ortho projection coordinates
float valueRemapping(float value, float low1, float high1, float low2, float high2)
{
	return low2 + (value - low1) * (high2 - low2) / (high1 - low1);
}	

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
	{
		cout << "right mouse button" << endl;
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		glm::vec2 result; //= valueRemapping(glm::vec2(xpos, ypos), 0, WINDOW_SIZE, 0, N);
		result.x = valueRemapping(xpos, 0, WINDOW_SIZE_WIDTH, 0, N);
		result.y = valueRemapping(ypos, 0, WINDOW_SIZE_HEIGHT, 0, N);
		// std::cout << xpos/double(WINDOW_SIZE) << ", " << 1.0 -ypos/double(WINDOW_SIZE)  << std::endl;
		std::cout << result.x << ", " << result.y << endl;
		dens_prev[IX(int(result.x), int((N+2) - result.y))] = 1.0;
		std::cout << "density sources value at " << result.x << ", "<< (N+2) - result.y << " : " << dens_prev[IX(int(result.x), int((N+2) - result.y))] << endl;
		// int index = 0;
		// for (int i=0;i<N+2;i++)
		// {
		// 	for (int j=0;j<N+2;j++)
		// 	{
		// 		gridPoints[index] = float(i);///float(N);
		// 		gridPoints[index+1] = float(j);///float(N);
		// 		gridPoints[index+2] = 0.0;  
		// 		gridPoints[index+3] += dens_prev[IX(i,j)]*20;
		// 		// std::cout << gridPoints[index]<< std::endl;
		// 		index = index + 4;
		// 	}
		// }
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		cout << "left mouse button" << endl;
	}
}

// not used yet
void drawGrid(float* gridValues, unsigned int VBO, unsigned int VAO, int shaderID)
{
	glBufferData(GL_ARRAY_BUFFER, sizeof(gridValues), gridValues, GL_STREAM_DRAW);
	
	GLint posAttrib = glGetAttribLocation(shaderID, "aPos");
	std::cout << posAttrib << std::endl;
	// iterpreting data from buffer 
	glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	// specify index of the vertex attribute to be enabled. --> 0
    glEnableVertexAttribArray(0);

	GLint densAttrib = glGetAttribLocation(shaderID, "densitySource");
	std::cout << densAttrib << std::endl;
	// iterpreting data from buffer 
	glVertexAttribPointer(densAttrib, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(3* sizeof(float)));
	// specify index of the vertex attribute to be enabled. --> 0
    glEnableVertexAttribArray(1);
}

// should this function is static
void allocate_data()
{
	// int size = (M+2)*(N+2)*(O+2);

	u			= (float *) malloc ( size*sizeof(float) );
	v			= (float *) malloc ( size*sizeof(float) );
	// w			= (float *) malloc ( size*sizeof(float) );
	u_prev		= (float *) malloc ( size*sizeof(float) );
	v_prev		= (float *) malloc ( size*sizeof(float) );
	// w_prev		= (float *) malloc ( size*sizeof(float) );
	dens		= (float *) malloc ( size*sizeof(float) );	
	dens_prev	= (float *) malloc ( size*sizeof(float) );


	if ( !dens || !dens_prev ) {
		fprintf ( stderr, "cannot allocate data\n" );
		// return ( 0 );
	}
	// initialise grid!
	for (int i=0;i<size;i++)
	{
		dens[i] = dens_prev[i] = u[i] = u_prev[i] = v[i] = v_prev[i] = 0.0f;
	}

	// return ( 1 );
}

// should this function be static ?
void free_data()
{
	for (int i=0 ; i<size ; i++ ) 
	{
		dens[i] = dens_prev[i] = 0.0f;
	}
}

int main()
{
	glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


	GLFWwindow* window = glfwCreateWindow(WINDOW_SIZE_WIDTH, WINDOW_SIZE_HEIGHT, "Fluid Sim", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	// tell GLFW to make the context of our window the main context on the current thread.
	glfwMakeContextCurrent(window);

	// window resize
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);  

	// initialize GLAD
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
    	std::cout << "Failed to initialize GLAD" << std::endl;
    	return -1;
	} 

	 // build and compile our shader program
    // ------------------------------------
    Shader ourShader("vertexShader.vs", "fragmentShader.fs"); 
	 
	// set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------

	allocate_data();

	// dens_prev = (float*) malloc((N+2)*(N+2)*sizeof(float));
	// for (int i=0;i<N+2;i++)
	// {
	// 	for (int j=0;j<N+2;j++)
	// 	{
	// 		dens_prev[IX(i,j)] = 0.0;
	// 	}
	// }


	// gridPoints = (float*) malloc((3*60*60)*sizeof(float));
	int index = 0;
	for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
			gridPoints[index] = float(i);///float(N);
			gridPoints[index+1] = float(j);///float(N);
			gridPoints[index+2] = 0.0;  
			gridPoints[index+3] = dens_prev[IX(i,j)];
			// std::cout << gridPoints[index]<< std::endl;
			index = index + 4;
		}
	}
	
	
	// float vertices[] = {
    //     -0.5f, -0.5f, 0.0f, // left  
    //     0.5f, -0.5f, 0.0f, // right 
    //     0.0f,  0.5f, 0.0f  // top   
    // }; 
    // unsigned int indices[] = {  // note that we start from 0!
    //     0, 1, 3,  // first Triangle
    //     1, 2, 3   // second Triangle
    // };
	// vel grid
	for (int i=0;i<size;i++)
	{
		u_prev[i] = v_prev[i] = 0.5; u[i] = v[i] = 0.5f;
	}


	int index1 = 0;
	for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
			vel_pos[index1] = float(i);
			vel_pos[index1+1] = float(j);
			vel_pos[index1+2] = 0.0f;
			// vel_pos[index1+3] = dens_prev[IX(i,j)];
			vel_pos[index1+3] = float(i) + u_prev[IX(i,j)];
			vel_pos[index1+4] = float(j) + v_prev[IX(i,j)];
			vel_pos[index1+5] = 0.0f;
			// cout << float(i) << ", " << float(j) << ", " << float(i) + 0.5f << ", " << float(j) + 0.5f << endl; 
			// vel_pos[index1+6] = dens_prev[IX(i,j)];	 
			
			index1 = index1 + 6;
		}
	}

	int index2 =0;
	for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
			dens_grid[index2] = float(i);
			dens_grid[index2+1] = float(j);
			dens_grid[index2+2] = 0.0f;
			dens_grid[index2+3] = dens_prev[IX(i,j)];
			index2 += 4;
		}
	}

	// generate buffers- vertex array object and vertex buffer object
	GLuint vaos[2];
    // unsigned int VBO, VAO, density_VBO;
	glGenVertexArrays(2, vaos);
	GLuint vbos[2];
    glGenBuffers(2, vbos);

	// Bind vertex array object first and then bind the vertex buffer objects
	// glBindVertexArray(VAO);

	// glBindBuffer(GL_ARRAY_BUFFER, VBO);
	
	// unbind
	// glBindBuffer(GL_ARRAY_BUFFER, 0);
	// glBindVertexArray(0);


	// glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	// render loop
	while(!glfwWindowShouldClose(window))
	{
		glBindVertexArray(vaos[0]);
		glBindBuffer(GL_ARRAY_BUFFER, vbos[0]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vel_pos), vel_pos, GL_STREAM_DRAW);

		GLint velAttrib = glGetAttribLocation(ourShader.ID, "aPos");
		// std::cout << posAttrib << std::endl;
		// iterpreting data from buffer 
		glVertexAttribPointer(velAttrib, 3, GL_FLOAT, GL_FALSE,  3 *sizeof(float), (void*)0);
		// specify index of the vertex attribute to be enabled. --> 0
		glEnableVertexAttribArray(0);

		glBindVertexArray(vaos[1]);
		glBindBuffer(GL_ARRAY_BUFFER, vbos[1]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(dens_grid), dens_grid, GL_STREAM_DRAW);

		GLint posAttrib = glGetAttribLocation(ourShader.ID, "aPos");
		// std::cout << posAttrib << std::endl;
		// iterpreting data from buffer 
		glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE,  4 *sizeof(float), (void*)0);
		// specify index of the vertex attribute to be enabled. --> 0
		glEnableVertexAttribArray(0);


		GLint densAttrib = glGetAttribLocation(ourShader.ID, "densitySource");
		// std::cout << posAttrib << std::endl;
		// iterpreting data from buffer 
		glVertexAttribPointer(densAttrib, 1, GL_FLOAT, GL_FALSE,  4 *sizeof(float), (void*)(3*sizeof(float)));
		// specify index of the vertex attribute to be enabled. --> 0
		glEnableVertexAttribArray(1);



		// GLint densAttrib = glGetAttribLocation(ourShader.ID, "densitySource");
		// std::cout << densAttrib << std::endl;
		// iterpreting data from buffer 
		// glVertexAttribPointer(densAttrib, 1, GL_FLOAT, GL_FALSE, 4* sizeof(float), (void*)(3* sizeof(float)));
		// specify index of the vertex attribute to be enabled. --> 0
		// glEnableVertexAttribArray(1);
	

		// begin simulation 
	

		// // input
		processInput(window);
		glfwSetMouseButtonCallback(window, mouse_button_callback);
		// rendering commands
		// ...
		// glClearColor is a state setting function, glClear is a state using function
		glClearColor(0.2f, 0.2f, 0.2f, 1.0f); 
		glClear(GL_COLOR_BUFFER_BIT);

		// glViewport(0,0,WINDOW_SIZE, WINDOW_SIZE);
		// this shader program/shaders are used for all rendering calls.
		// // glUseProgram(shaderProgram);
		ourShader.use();

		glm::mat4 trans = glm::mat4(1.0f);
		trans = glm::translate(trans, glm::vec3(-0.5f, -0.5f, 0.0f));
		unsigned int transformMatrixLocation = glGetUniformLocation(ourShader.ID, "transform");
		// arg 2- how many matrices, arg3-> should matrix be transposed
		glUniformMatrix4fv(transformMatrixLocation, 1, GL_FALSE, glm::value_ptr(trans));

		glm::mat4 projection = glm::ortho(-10.0f, 110.0f, -1.0f, 110.0f, -1.0f, 100.0f);

		unsigned int projectionMatrixLocation = glGetUniformLocation(ourShader.ID, "projection");
		// arg 2- how many matrices, arg3-> should matrix be transposed
		glUniformMatrix4fv(projectionMatrixLocation, 1, GL_FALSE, glm::value_ptr(projection));

		glBindVertexArray(vaos[0]); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
		// glPointSize(4.0f);
		glLineWidth(1.0f);
	    glDrawArrays(GL_LINES,  0, 2*size);
		// glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		// glDrawElements(GL_LINES, 10, GL_UNSIGNED_INT, 0);
		// check and call events and swap buffers

		glBindVertexArray(vaos[1]); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
		glPointSize(4.0f);
		// glLineWidth(1.0f);
	    glDrawArrays(GL_POINTS,  0, size);
		// glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
		// glDrawElements(GL_LINES, 10, GL_UNSIGNED_INT, 0);
		// check and call events and swap buffers

		glfwSwapBuffers(window);
		glfwPollEvents(); 
		// drawGrid(gridPoints, VBO, VAO,  ourShader.ID);

		dens_step (N, dens, dens_prev, u, v, diff, dt , vel_pos, dens_grid);
	}

	// optional: de-allocate all resources once they've outlived their purpose
    glDeleteVertexArrays(2, vaos);
    glDeleteBuffers(2, vbos);

	// glfw: terminate, clearing all previously allocated GLFW resources.
	glfwTerminate();
	return 0;
}
