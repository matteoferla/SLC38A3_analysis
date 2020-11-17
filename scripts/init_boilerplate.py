def make_option_string(**options):
    """
    This just converts the key:value pairs to a command line string for the pyrosetta init.

    >> make_option_string(no_optH=False, ex1=None, remodel=dict(blueprint='mod.blue'))

    Bools are converted,
    None results in a value argument,
    Dictionaries are converted to xx:xx type arguments

    Also... Full option list: https://www.rosettacommons.org/docs/latest/full-options-list
    """

    def format_inner(v):
        if v is None:
            return ''
        elif v in (False, True):
            return ' ' + str(v).lower()
        else:
            return ' ' + str(v)

    def format_outer(k):
        if isinstance(options[k], dict):
            return ' '.join([f'-{k}:{k2}' + format_inner(v2) for k2, v2 in options[k].items()])
        else:
            return '-' + k + format_inner(options[k])

    args = [format_outer(k) for k in options]
    return ' '.join(args)

