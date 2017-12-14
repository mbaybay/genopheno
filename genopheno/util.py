import time
import os


def expand_path(path):
    """
    Expands a file path. This handles users and environmental variables.
    :param path: The path the expand
    :return: The expanded path
    """
    new_path = os.path.expanduser(path)
    return os.path.expandvars(new_path)


def clean_output(output_dir):
    """
    Creates the output directory if it does not exist and removes old files if it does.
    :param output_dir:
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        for f in os.listdir(output_dir):
            os.remove(os.path.join(output_dir, f))


def timed_invoke(action, method):
    """
    Invokes a method. Prints the start and finish with the invocation time.
    :param action: The string describing the method action
    :param method: The method to invoke
    :return: The return of the method
    """
    print 'Started {}...'.format(action)
    start = time.time()
    try:
        output = method()
        print 'Finished {} in {} seconds'.format(action, int(time.time() - start))
        return output
    except Exception:
        print 'Exception while {} after {} seconds'.format(action, int(time.time() - start))
        raise
