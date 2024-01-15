import sys
sys.path.append('../../pebble_factory/')
import triso_pebble
import openmc 
import numpy as np
import openmc.model


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
    bounds = [lowerLeft[0],lowerLeft[1],lowerLeft[2],
        upperRight[0],upperRight[1],upperRight[2]]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:],
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
    triso_pebble_keff_test()

