/**
 * common.glsl
 * Common types and functions used for ray tracing.
 */

const float pi = 3.14159265358979;
const float epsilon = 0.001;

struct Ray {
    vec3 o;     // origin
    vec3 d;     // direction - always set with normalized vector
    float t;    // time, for motion blur
};

Ray createRay(vec3 o, vec3 d, float t)
{
    Ray r;
    r.o = o;
    r.d = d;
    r.t = t;
    return r;
}

Ray createRay(vec3 o, vec3 d)
{
    return createRay(o, d, 0.0);
}

vec3 pointOnRay(Ray r, float t)
{
    return r.o + r.d * t;
}

float gSeed = 0.0;

uint baseHash(uvec2 p)
{
    p = 1103515245U * ((p >> 1U) ^ (p.yx));
    uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
    return h32 ^ (h32 >> 16);
}

float hash1(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    return float(n) / float(0xffffffffU);
}

vec2 hash2(inout float seed) {
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1,seed += 0.1)));
    uvec2 rz = uvec2(n, n * 48271U);
    return vec2(rz.xy & uvec2(0x7fffffffU)) / float(0x7fffffff);
}

vec3 hash3(inout float seed)
{
    uint n = baseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
    uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
    return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

float rand(vec2 v)
{
    return fract(sin(dot(v.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 toLinear(vec3 c)
{
    return pow(c, vec3(2.2));
}

vec3 toGamma(vec3 c)
{
    return pow(c, vec3(1.0 / 2.2));
}

vec2 randomInUnitDisk(inout float seed) {
    vec2 h = hash2(seed) * vec2(1.0, 6.28318530718);
    float phi = h.y;
    float r = sqrt(h.x);
	return r * vec2(sin(phi), cos(phi));
}

vec3 randomInUnitSphere(inout float seed)
{
    vec3 h = hash3(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
    float phi = h.y;
    float r = pow(h.z, 1.0/3.0);
	return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

vec3 randomUnitVector(inout float seed) //to be used in diffuse reflections with distribution cosine
{
    return(normalize(randomInUnitSphere(seed)));
}

struct Camera
{
    vec3 eye;
    vec3 u, v, n;
    float width, height;
    float lensRadius;
    float planeDist, focusDist;
    float time0, time1;
};

Camera createCamera(
    vec3 eye,
    vec3 at,
    vec3 worldUp,
    float fovy,
    float aspect,
    float aperture,  //diametro em multiplos do pixel size
    float focusDist,  //focal ratio
    float time0,
    float time1)
{
    Camera cam;
    if(aperture == 0.0) cam.focusDist = 1.0; //pinhole camera then focus in on vis plane
    else cam.focusDist = focusDist;
    vec3 w = eye - at;
    cam.planeDist = length(w);
    cam.height = 2.0 * cam.planeDist * tan(fovy * pi / 180.0 * 0.5);
    cam.width = aspect * cam.height;

    cam.lensRadius = aperture * 0.5 * cam.width / iResolution.x;  //aperture ratio * pixel size; (1 pixel=lente raio 0.5)
    cam.eye = eye;
    cam.n = normalize(w);
    cam.u = normalize(cross(worldUp, cam.n));
    cam.v = cross(cam.n, cam.u);
    cam.time0 = time0;
    cam.time1 = time1;
    return cam;
}

Ray getRay(Camera cam, vec2 pixel_sample)  //rnd pixel_sample viewport coordinates
{
    vec2 ls = cam.lensRadius * randomInUnitDisk(gSeed);  //ls - lens sample for DOF
    float time = cam.time0 + hash1(gSeed) * (cam.time1 - cam.time0);
    
	vec3 ray_dir = cam.u * cam.width * ((pixel_sample.x) / iResolution.x - 0.5) +
		cam.v * cam.height * ((pixel_sample.y) / iResolution.y - 0.5) -
		cam.n *(cam.planeDist);
    
    return createRay(cam.eye, normalize(ray_dir), time);
}

// MT_ material type
#define MT_DIFFUSE 0
#define MT_METAL 1
#define MT_DIALECTRIC 2

struct Material
{
    int type;
    vec3 albedo;  //diffuse color
    vec3 specColor;  //the color tint for specular reflections. for metals and opaque dieletrics like coloured glossy plastic
    vec3 emissive; //
    float roughness; // controls roughness for metals. It can be used for rough refractions
    float refIdx; // index of refraction for dialectric
    vec3 refractColor; // absorption for beer's law
};

Material createDiffuseMaterial(vec3 albedo)
{
    Material m;
    m.type = MT_DIFFUSE;
    m.albedo = albedo;
    m.specColor = vec3(0.0);
    m.roughness = 1.0;  //ser usado na iluminação direta
    m.refIdx = 1.0;
    m.refractColor = vec3(0.0);
    m.emissive = vec3(0.0);
    return m;
}

Material createMetalMaterial(vec3 specClr, float roughness)
{
    Material m;
    m.type = MT_METAL;
    m.albedo = vec3(0.0);
    m.specColor = specClr;
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

Material createDialectricMaterial(vec3 refractClr, float refIdx, float roughness)
{
    Material m;
    m.type = MT_DIALECTRIC;
    m.albedo = vec3(0.0);
    m.specColor = vec3(0.04);
    m.refIdx = refIdx;
    m.refractColor = refractClr;  
    m.roughness = roughness;
    m.emissive = vec3(0.0);
    return m;
}

struct HitRecord
{
    vec3 pos;
    vec3 normal;
    float t;            // ray parameter
    Material material;
};

bool customRefract(vec3 v, vec3 n, float niOverNt, out vec3 refracted) {
  vec3 uv = normalize(v);
  float dt = dot(uv, n);
  float discriminant = 1.0 - niOverNt * niOverNt * (1.0 - dt * dt);
  if (discriminant > 0.0) {
    refracted = niOverNt * (uv - n * dt) - n * sqrt(discriminant);
    return true; // Refraction occurred
  } else {
    refracted = vec3(0.0);
    return false; // Total internal reflection
  }
}

vec3 metalSchlick(float cosTheta, vec3 specColor) {
  return specColor + (vec3(1.0) - specColor) * pow(1.0 - cosTheta, 5.0);
}

float schlick(float cosine, float refIdx)
{
    float r0 = (1.0 - refIdx) / (1.0 + refIdx);
    r0 = r0 * r0;
   
	//return r0 + (1.0 - r0) * pow(1.0 - cosine, 5.0);
    return r0 + (1.0 - r0) * pow(clamp(1.0 - cosine, 0.0, 1.0), 5.0);
}

bool scatter(Ray rIn, HitRecord rec, out vec3 atten, out Ray rScattered)
{
    vec3 recPos;
    if(rec.material.type == MT_DIFFUSE)
    {
        //INSERT CODE HERE, arranjar para não ficar igual ao prof
        vec3 normal = rec.normal;
        if(dot(rIn.d, rec.normal) > 0.0)
            normal *= -1.0;
        vec3 target = rec.pos + normal + randomUnitVector(gSeed);
        rScattered = createRay(rec.pos + epsilon * normal, normalize(target - rec.pos), rIn.t);
        //float cosTheta = dot( n, l);
        atten = rec.material.albedo * max(dot(rScattered.d, rec.normal), 0.0) / pi;
        return true;
    }
    if(rec.material.type == MT_METAL)
    {
       //INSERT CODE HERE, consider fuzzy reflections
        vec3 refl = reflect(rIn.d, rec.normal);
        recPos = epsilon + rec.pos * rec.normal;
        vec2 perturbation = rec.material.roughness * randomInUnitDisk(gSeed);
        vec3 scatteredDirection = normalize(refl + vec3(perturbation, 0.0));
        rScattered= createRay(recPos, scatteredDirection, rIn.t);
        atten = metalSchlick(-dot(rIn.d, rec.normal), rec.material.specColor);
        //atten = rec.material.specColor;
        return true;
    }
    if(rec.material.type == MT_DIALECTRIC)
    {
        atten = vec3(1.0);
        vec3 outwardNormal;
        float niOverNt;
        float cosine;

        cosine = -dot(rIn.d, rec.normal);

        if(dot(rIn.d, rec.normal) > 0.0) //hit inside
        {
            outwardNormal = -rec.normal;
            niOverNt = rec.material.refIdx;
            atten = exp(-rec.material.refractColor*rec.t); //beer law
            //cosine = refraction cosine for schlick; 
            //atten = apply Beer's law by using rec.material.refractColor
        }
        else  //hit from outside
        {
            outwardNormal = rec.normal;
            niOverNt = 1.0 / rec.material.refIdx;
           // cosine = -dot(rIn.d, rec.normal); 
        }

        //Use probabilistic math to decide if scatter a reflected ray or a refracted ray
        vec3 refracted;
        float reflectProb;

        if(customRefract(rIn.d, outwardNormal, niOverNt, refracted)) //no total internal reflection
        {
            if(niOverNt > 1.0) {
                cosine = sqrt(1.0-niOverNt*niOverNt*(1.0 * cosine*cosine));
            }
            reflectProb = schlick(cosine, rec.material.refIdx);
        }
        else
        {
            reflectProb = 1.0;
        }
        //if no total reflection  reflectProb = schlick(cosine, rec.material.refIdx);  
        //else reflectProb = 1.0;

        if( hash1(gSeed) < reflectProb)  //Reflection
            rScattered = createRay(recPos, reflect(rIn.d, rec.normal));
          // atten *= vec3(reflectProb); not necessary since we are only scattering reflectProb rays and not all reflected rays
        
        else  //Refraction
            rScattered = createRay(recPos, refracted);
           // atten *= vec3(1.0 - reflectProb); not necessary since we are only scattering 1-reflectProb rays and not all refracted rays
    
        return true;
    }
    return false;
}

struct Triangle {vec3 a; vec3 b; vec3 c; };

Triangle createTriangle(vec3 v0, vec3 v1, vec3 v2)
{
    Triangle t;
    t.a = v0; t.b = v1; t.c = v2;
    return t;
}

bool hit_triangle(Triangle triangle, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //calculate a valid t and normal
    vec3 normal;
    float t;
    vec3 edge1, edge2, edge3, tvec, pvec, qvec;

	float det, inv_det;

	edge1 = triangle.b - triangle.a;
	edge2 = triangle.c - triangle.a;

    normal = cross(edge1, edge2);
	normalize(normal);

	pvec = cross(r.d, edge2);

	det = dot(edge1, pvec);

	if (det < epsilon) return false;

	tvec = r.o - triangle.a;

	float u = dot(tvec, pvec);
	if (u < 0.0 || u > det) return false;

	qvec = cross(tvec, edge1);

	float v = dot(r.d, qvec);
	if (v < 0.0 || u + v > det) return false;

	t = dot(edge2, qvec);
	inv_det = 1.f / det;
	t *= inv_det;
	u *= inv_det;
	v *= inv_det;
  
    if(t < tmax && t > tmin)
    {
        rec.t = t;
        rec.normal = normal;
        rec.pos = pointOnRay(r, rec.t);
        return true;
    }
    return false;
}


struct Sphere
{
    vec3 center;
    float radius;
};

Sphere createSphere(vec3 center, float radius)
{
    Sphere s;
    s.center = center;
    s.radius = radius;
    return s;
}


struct MovingSphere
{
    vec3 center0, center1;
    float radius;
    float time0, time1;
};

MovingSphere createMovingSphere(vec3 center0, vec3 center1, float radius, float time0, float time1)
{
    MovingSphere s;
    s.center0 = center0;
    s.center1 = center1;
    s.radius = radius;
    s.time0 = time0;
    s.time1 = time1;
    return s;
}

vec3 center(MovingSphere mvsphere, float time)
{
    //return mvsphere.center0 + ((time - mvsphere.time0) / (mvsphere.time1 - mvsphere.time0)* center1-center0?);
    return vec3(0.0,0.0,0.0);
}


/*
 * The function naming convention changes with these functions to show that they implement a sort of interface for
 * the book's notion of "hittable". E.g. hit_<type>.
 */

bool hit_sphere(Sphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    //calculate a valid t and normal
    float t;
    vec3 normal;
    vec3 oc = s.center - r.o;
	float b = dot(r.d, oc);
	float c = dot(oc, oc) - pow(s.radius, 2.0);

    if (c < 0.0) {
		t = b + sqrt(pow(b, 2.0) - c);
	}
	else {
		if (b <= 0.0) return false;
		else {
			if ((pow(b, 2.0) - c) <= 0.0) return false;
			else {
				t = b - sqrt(pow(b, 2.0) - c);
			}
		}
	}

    normal = rec.pos - s.center;
	
    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);
        rec.normal = normal;
        return true;
    }
    else return false;
}

bool hit_movingSphere(MovingSphere s, Ray r, float tmin, float tmax, out HitRecord rec)
{
    float B, C, delta;
    bool outside;
    float t;
    vec3 normal;

   // vec3 sphereCen;

     //INSERT YOUR CODE HERE
     //Calculate the moving center
    //calculate a valid t and normal
	
    if(t < tmax && t > tmin) {
        rec.t = t;
        rec.pos = pointOnRay(r, rec.t);
        rec.normal = normal;
        return true;
    }
    else return false;
}

struct pointLight {
    vec3 pos;
    vec3 color;
};

pointLight createPointLight(vec3 pos, vec3 color) 
{
    pointLight l;
    l.pos = pos;
    l.color = color;
    return l;
}
