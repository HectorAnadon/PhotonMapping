/*********************************************************************************
Copyright (C) 2014 Adrian Jarabo (ajarabo@unizar.es)
Copyright (C) 2014 Diego Gutierrez (diegog@unizar.es)
All rights reserved.

This is an educational Ray Tracer developed for the course 'Informatica Grafica'
(Computer Graphics) tought at Universidad de Zaragoza (Spain). As such, it does not
intend to be fast or general, but just to provide an educational tool for undergraduate
students.

This software is provided as is, and any express or implied warranties are disclaimed.
In no event shall copyright holders be liable for any damage.
**********************************************************************************/
#include "PhotonMapping.h"
#include "World.h"
#include "Intersection.h"
#include "Ray.h"
#include "BSDF.h"
#include "math.h"
#include <iostream>

//*********************************************************************
// Compute the photons by tracing the Ray 'r' from the light source
// through the scene, and by storing the intersections with matter
// in the lists 'xx_photons', storing the diffuse (global) and caustic
// photons respectively. For efficiency, both are computed at the same
// time, since computing them separately would result into a lost of
// several samples marked as caustic or diffuse.
// Same goes with the boolean 'direct', that specifies if direct 
// photons (from light to surface) are being stored or not. 
// The initial traced photon has energy defined by the tristimulus
// 'p', that accounts for the emitted power of the light source.
// The function will return true when there are more photons (caustic
// or diffuse) to be shot, and false otherwise.
//---------------------------------------------------------------------
bool PhotonMapping::trace_ray(const Ray& r, const Vector3 &p,
	std::list<Photon> &global_photons, std::list<Photon> &caustic_photons, bool direct)
{

	//Check if max number of shots done...
	if (++m_nb_current_shots > m_max_nb_shots)
	{
		return false;
	}

	// Compute irradiance photon's energy
	Vector3 energy(p);

	Ray photon_ray(r);
	photon_ray.shift();

	bool is_caustic_particle = false;

	//Iterate the path
	while (1)
	{
		// Throw ray and update current_it
		Intersection it;
		world->first_intersection(photon_ray, it);

		if (!it.did_hit())
			break;

		//Check if has hit a delta material...
		if (it.intersected()->material()->is_delta())
		{
			// If delta material, then is caustic...
			// Don't store the photon!
			is_caustic_particle = true;
		}
		else if (photon_ray.get_level() > 0 || direct)
		{
			//If non-delta material, store the photon!
			if (is_caustic_particle)
			{
				//If caustic particle, store in caustics
				if (caustic_photons.size() < m_nb_caustic_photons)
					caustic_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy));
			}
			else
			{
				//If non-caustic particle, store in global
				if (global_photons.size() < m_nb_global_photons)
					global_photons.push_back(Photon(it.get_position(), photon_ray.get_direction(), energy));
			}
			is_caustic_particle = false;
		}

		Real pdf;

		Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
		Real avg_surf_albedo = surf_albedo.avg();

		Real epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);
		while (epsilon2 < 0.)
			epsilon2 = static_cast<Real>(rand()) / static_cast<Real>(RAND_MAX);

		if (epsilon2 > avg_surf_albedo || photon_ray.get_level() > 20)
			break;

		// Random walk's next step
		// Get sampled direction plus pdf, and update attenuation
		it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf);

		// Shade...
		energy = energy*surf_albedo;
		if (!it.intersected()->material()->is_delta())
			energy *= dot_abs(it.get_normal(), photon_ray.get_direction()) / 3.14159;

		energy = energy / (pdf*avg_surf_albedo);
	}

	if (caustic_photons.size() == m_nb_caustic_photons &&
		global_photons.size() == m_nb_global_photons)
	{
		m_max_nb_shots = m_nb_current_shots - 1;
		return false;
	}

	return true;
}

double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f* (fMax - fMin);
}

/*Return a random Vector3d*/
Vector3 randomDirection() {
	double o(fRand(0, 2 * 3.141593));
	double t(fRand(0, 2 * 3.141593));
	double x = cos(t) * sin(o);
	double y = cos(t) * cos(o);
	double z = sin(t);
	Vector3 direction(x, y, z);
	return direction;
}

//*********************************************************************
// TODO: Implement the preprocess step of photon mapping,
// where the photons are traced through the scene. To do it,
// you need to follow these steps for each shoot:
//  1 - Sample a world's light source in the scene to create
//		the initial direct photon from the light source.
//	2 - Trace the photon through the scene storing the inter-
//		sections between the photons and matter. You can use
//		the function 'trace_ray' for this purpose.
//	3 - Finally, once all the photons have been shot, you'll
//		need to build the photon maps, that will be used later
//		for rendering. 
//		NOTE: Careful with function
//---------------------------------------------------------------------
void PhotonMapping::preprocess() {
	std::list<Photon> global;
	std::list<Photon> caustic;
	Ray* r;

	//Loop for all lights
	for (int i = 0; i < world->nb_lights(); i++) {
		Vector3 light_position = world->light(i).get_position();
		Vector3 light_intensity = world->light(i).get_intensities();
		LightSource* light_source = new PointLightSource(world, light_position, light_intensity);
		do
		{
			//Get random direction to shoot a photon
			Vector3 direction = randomDirection();
			r = new Ray(light_position, direction);
			//Until all photons are done (20 rebounds per photon or not hit/glass/mirror)
		} while (trace_ray(*r, light_intensity, global, caustic, false));
	}

	//Insert global_photons list into kdtree
	while (global.size() > 0) {
		Photon photon = global.front();
		global.pop_front();
		m_global_map.store(std::vector<Real>(photon.position.data, photon.position.data + 3), photon);
	}
	m_global_map.balance();

	//Insert caustic_photons list into kdtree
	if (caustic.size() > 0) {
		while (caustic.size() > 0) {
			Photon photon = caustic.front();
			caustic.pop_front();
			m_caustics_map.store(std::vector<Real>(photon.position.data, photon.position.data + 3), photon);
		}
		m_caustics_map.balance();
	}

}

//*********************************************************************
// TODO: Implement the function that computes the rendering equation 
// using radiance estimation with photon mapping, using the photon
// maps computed as a proprocess. Note that you will need to handle
// both direct and global illumination, together with recursive the 
// recursive evaluation of delta materials. For an optimal implemen-
// tation you should be able to do it iteratively.
// In principle, the class is prepared to perform radiance estimation
// using k-nearest neighbors ('m_nb_photons') to define the bandwidth
// of the kernel.
//---------------------------------------------------------------------
Vector3 PhotonMapping::shade(Intersection &it0)const {
	Intersection it(it0);

	int current_bounds = 0;
	Ray transmitted_ray;
	while (it.intersected()->material()->is_delta() && current_bounds < maximum_bounds) {
		Real pdf;	//Useless
		it.intersected()->material()->get_outgoing_sample_ray(it, transmitted_ray, pdf);
		transmitted_ray.shift();
		world->first_intersection(transmitted_ray, it);
		current_bounds++;
	}
	Vector3 intersection = it.get_position();
	Vector3 normal = it.get_normal();

	Vector3 L(world->get_background());

	//Ambiental light
	//Vector3 L = it.intersected()->material()->get_albedo(it) * world->get_ambient();

	//Direct light
	//For each light
	for (int i = 0; i < world->nb_lights(); i++) {
		Vector3 light_position = world->light(i).get_position();
		Vector3 light_intensity = world->light(i).get_intensities();
		LightSource* light_source = new PointLightSource(world, light_position, light_intensity);

		Vector3 r_light = Vector3(0) - light_source->get_incoming_direction(intersection);
		//If direct light
		if (light_source->is_visible(intersection)) {
			Vector3 id = light_source->get_incoming_light(intersection);
			Vector3 kd = it.intersected()->material()->get_albedo(it);
			Vector3 r_light_norm = r_light.normalize();
			Real cos = r_light_norm.dot(normal);
			L += kd * id * cos;		//Diffuse

			/*Vector3 is = light_source->get_incoming_light(intersection);
			Vector3 ks = it.intersected()->material()->get_albedo(it);
			Vector3 R = r_light.reflect(normal);
			Vector3 V = it.get_ray().get_direction();
			R = R.normalize();
			V = V.normalize();
			cos = R.dot(V);
			if (cos > 0.0) {
				L += ks * is * pow(cos, 100);	//Specular
			}*/
		}
	}


	std::vector<const KDTree<Photon, 3>::Node*> global_photons;
	Real max_distance = 100;
	m_global_map.find(std::vector<Real>(it.get_position().data, it.get_position().data + 3), m_nb_photons, global_photons, max_distance);

	//cout << global_photons.size();

	std::vector<const KDTree<Photon, 3>::Node*> causics_photons;
	m_caustics_map.find(std::vector<Real>(it.get_position().data, it.get_position().data + 3), m_nb_photons, causics_photons, max_distance);


	if (global_photons.size() > 1) {
		//cout << "Global " << global_photons.size() << "\n";

		Vector3 color_bleeding(0);
		Real max_radious = 0.0;

		for (long i = 0; i<global_photons.size(); ++i) {
			KDTree<Photon, 3>::Node node = *global_photons[i];
			Photon photon = node.data();
			Vector3 photon_position = photon.position;

			Real photon_distance = sqrt(pow(photon_position[0] - intersection[0], 2) + pow(photon_position[1] - intersection[1], 2) + 
				pow(photon_position[2] - intersection[2], 2));
			if (max_radious < photon_distance) {
				max_radious = photon_distance;
			}

			color_bleeding += photon.flux * it.intersected()->material()->get_albedo(it);

		}

		// Circle area
		Real area = 3.141593 * pow(max_radious, 2);
		//L += color_bleeding / (area*global_photons.size());
		L += color_bleeding / (area*m_nb_global_photons);
		//L += color_bleeding / (area);
	}


	if (causics_photons.size() > 1) {
		//cout << "Caustic " << causics_photons.size() << "\n";

		Vector3 caustic(0);
		Real max_radious = 0.0;

		for (long i = 0; i<causics_photons.size(); ++i) {
			KDTree<Photon, 3>::Node node = *causics_photons[i];
			Photon photon = node.data();
			Vector3 photon_position = photon.position;

			Real photon_distance = sqrt(pow(photon_position[0] - intersection[0], 2) + pow(photon_position[1] - intersection[1], 2) + 
				pow(photon_position[2] - intersection[2], 2));
			if (max_radious < photon_distance) {
				max_radious = photon_distance;
			}

			caustic += photon.flux * it.intersected()->material()->get_albedo(it);

			//cout << "Photon " << photon.direction.data[0] << "\n";
		}

		Real area = 3.141593 * pow(max_radious, 2);
		//L += caustic / (area*causics_photons.size());
		L += caustic / (area*m_nb_caustic_photons);
		//L += caustic / (area);

	}

	//cout << m_nb_photons << "\n";
	return L;

	//**********************************************************************
	// The following piece of code is included here for two reasons: first
	// it works as a 'hello world' code to check that everthing compiles 
	// just fine, and second, to illustrate some of the functions that you 
	// will need when doing the work. Goes without saying: remove the 
	// pieces of code that you won't be using.
	//
	unsigned int debug_mode = 1;

	switch (debug_mode)
	{
	case 1:
		// ----------------------------------------------------------------
		// Display Albedo Only
		L = it.intersected()->material()->get_albedo(it);
		break;
	case 2:
		// ----------------------------------------------------------------
		// Display Normal Buffer
		L = it.get_normal();
		break;
	case 3:
		// ----------------------------------------------------------------
		// Display whether the material is specular (or refractive) 
		L = Vector3(it.intersected()->material()->is_delta());
		break;

	case 4:
		// ----------------------------------------------------------------
		// Display incoming illumination from light(0)
		L = world->light(0).get_incoming_light(it.get_position());
		break;

	case 5:
		// ----------------------------------------------------------------
		// Display incoming direction from light(0)
		L = world->light(0).get_incoming_direction(it.get_position());
		break;

	case 6:
		// ----------------------------------------------------------------
		// Check Visibility from light(0)
		if (world->light(0).is_visible(it.get_position()))
			L = Vector3(1.);
		break;
	}
	// End of exampled code
	//**********************************************************************

	return L;
}