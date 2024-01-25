import sys
sys.path.append('../../pebble_factory/')
import openmc.deplete


def postprocess_depletion():

    results = openmc.deplete.Results("./depletion_results.h5")
    time, k = results.get_keff()

    time /= (24 * 60 * 60)  # convert back to days from seconds

    from matplotlib import pyplot

    pyplot.errorbar(time, k[:, 0], yerr=k[:, 1])
    pyplot.xlabel("Time [d]")
    pyplot.ylabel("$k_{eff}\pm \sigma$");
    pyplot.savefig("depletion_results.jpg")



if __name__ == '__main__':

    import doctest
    doctest.testmod(verbose=True)

    # triso_pebble_keff_test()
    postprocess_depletion()

