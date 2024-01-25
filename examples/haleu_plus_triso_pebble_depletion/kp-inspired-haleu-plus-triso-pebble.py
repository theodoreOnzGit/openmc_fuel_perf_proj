import sys
sys.path.append('../../pebble_factory/')
import triso_pebble
import openmc
import numpy as np
import openmc.model
import openmc.deplete


def triso_pebble_depletion_test(light_uranium_u236_fraction):
    # import triso particle constructor
    from triso_haleu_plus import TrisoParticlesFactory
    triso_factory = TrisoParticlesFactory()
    fuel_temp_deg_kelvin = 600.0
    coolant_temp_deg_kelvin = 600.0
    pebble_shell_thickness_cm = 0.1
    flibe_lithium_7_enrichment_atom_percent = 99.995
    light_uranium_enrichment_atom_percent = 19.9

    # create materials for fresh triso pebble with normal haleu
    pebble_materials = triso_factory.build_haleu_plus_fhr_materials(
            light_uranium_enrichment_atom_percent,
            light_uranium_u236_fraction,
            flibe_lithium_7_enrichment_atom_percent,
            coolant_temp_deg_kelvin)

    # create triso pebble universe (fill)
    from triso_pebble import TrisoPebbleFactory
    pebble_factory = TrisoPebbleFactory()

    pebble_univ = pebble_factory.buildNonAnnularPebbleUniv(
            pebble_shell_thickness_cm,
            pebble_materials,
            fuel_temp_deg_kelvin,
            coolant_temp_deg_kelvin)

    # construct one single pebble
    test_sphere = openmc.Sphere(r=3.0)
    test_sphere.boundary_type = 'reflective'
    fuel_region = -test_sphere

    test_pebble = openmc.Cell()
    test_pebble.region = fuel_region
    test_pebble.fill = pebble_univ
    root_universe = openmc.Universe(cells=[test_pebble])
    geometry = openmc.Geometry(root_universe)

    # make geometry and materials
    geometry = openmc.Geometry(root_universe)
    materials = openmc.Materials()

    # for fuel, I need to specify its volume
    # note that this is only a placeholder,
    # i assume 20,000 triso particles with radius of 0.0215 cm
    fuel_vol_single_triso = (
            4.0/3.0 * np.float_power(0.0215, 3.0)
            * np.pi)

    fuel_vol_total_triso_20k = 20000 * fuel_vol_single_triso
    pebble_materials.fuel.volume = fuel_vol_total_triso_20k
    materials.append(pebble_materials.fuel)
    materials.append(pebble_materials.buff)
    materials.append(pebble_materials.PyC1)
    materials.append(pebble_materials.PyC2)
    materials.append(pebble_materials.SiC)
    materials.append(pebble_materials.graphite)
    materials.append(pebble_materials.flibe)

    # settings

    settings = openmc.Settings()
    settings.batches = 200
    settings.inactive = 10
    settings.particles = 5000
    lowerLeft, upperRight = pebble_univ.bounding_box
    bounds = [lowerLeft[0], lowerLeft[1], lowerLeft[2],
              upperRight[0], upperRight[1], upperRight[2]]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],
                                    only_fissionable=True)
    settings.source = openmc.Source(space=uniform_dist)

    # shannon entropy mesh
    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left, entropy_mesh.upper_right = geometry.bounding_box
    entropy_mesh.dimension = (8, 8, 8)
    settings.entropy_mesh = entropy_mesh

    # import depletion chain
    # uses a pwr neutron spectrum I suppose it approximates
    # a FHR better than using a SFR spectrum
    # probably need to generate my own spectrum in future
    chain = openmc.deplete.Chain.from_xml("../chain_endfb71_pwr.xml")

    # construct model
    model = openmc.Model(geometry=geometry, settings=settings)
    operator = openmc.deplete.CoupledOperator(
            model, "../chain_endfb71_pwr.xml")

    # power density of 174 W/cm, same as in the pin cell analysis tutorial 
    power_density_watts_per_cm = 174

    # time steps
    time_step_days = [30] * 6

    # use Euler integration scheme
    integrator = openmc.deplete.PredictorIntegrator(
            operator,
            time_step_days,
            power_density_watts_per_cm,
            timestep_units='d')

    # start depletion
    integrator.integrate()

def triso_pebble_keff_test():
    # import triso particle constructor
    from triso import TrisoParticlesFactory
    triso_factory = TrisoParticlesFactory()
    fuel_temp_deg_kelvin = 600.0
    coolant_temp_deg_kelvin = 600.0
    pebble_shell_thickness_cm = 0.1
    flibe_lithium_7_enrichment_atom_percent = 99.995
    uranium_235_enrichment_atom_percent = 19.9

    # create materials for fresh triso pebble with normal haleu
    pebble_materials = triso_factory.build_fhr_materials(
            uranium_235_enrichment_atom_percent,
            flibe_lithium_7_enrichment_atom_percent,
            coolant_temp_deg_kelvin)

    # create triso pebble universe (fill)
    from triso_pebble import TrisoPebbleFactory
    pebble_factory = TrisoPebbleFactory()

    pebble_univ = pebble_factory.buildNonAnnularPebbleUniv(
            pebble_shell_thickness_cm,
            pebble_materials,
            fuel_temp_deg_kelvin,
            coolant_temp_deg_kelvin)

    # construct one single pebble
    test_sphere = openmc.Sphere(r=3.0)
    test_sphere.boundary_type = 'reflective'
    fuel_region = -test_sphere

    test_pebble = openmc.Cell()
    test_pebble.region = fuel_region
    test_pebble.fill = pebble_univ
    root_universe = openmc.Universe(cells=[test_pebble])
    geometry = openmc.Geometry(root_universe)

    # make geometry and materials
    geometry = openmc.Geometry(root_universe)
    materials = openmc.Materials()
    materials.append(pebble_materials.fuel)
    materials.append(pebble_materials.buff)
    materials.append(pebble_materials.PyC1)
    materials.append(pebble_materials.PyC2)
    materials.append(pebble_materials.SiC)
    materials.append(pebble_materials.graphite)
    materials.append(pebble_materials.flibe)

    # settings

    settings = openmc.Settings()
    settings.batches = 200
    settings.inactive = 50
    settings.particles = 5000
    lowerLeft, upperRight = pebble_univ.bounding_box
    bounds = [lowerLeft[0], lowerLeft[1], lowerLeft[2],
              upperRight[0], upperRight[1], upperRight[2]]
    uniform_dist = openmc.stats.Box(
            bounds[:3],
            bounds[3:],
            only_fissionable=True)
    settings.source = openmc.Source(space=uniform_dist)

    # shannon entropy mesh
    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left, entropy_mesh.upper_right = geometry.bounding_box
    entropy_mesh.dimension = (8, 8, 8)
    settings.entropy_mesh = entropy_mesh

    # export to xml
    materials.export_to_xml()
    geometry.export_to_xml()
    settings.export_to_xml()

    # run openmc
    openmc.run()


if __name__ == '__main__':

    import doctest
    doctest.testmod(verbose=True)

    # triso_pebble_keff_test()
    triso_pebble_depletion_test(light_uranium_u236_fraction=0.5)
