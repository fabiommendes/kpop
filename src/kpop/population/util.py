def id_label_from_parents(l1, l2):
    """
    Creates a new id label from parents id labels.
    """

    if l1 == l2:
        return l1
    elif l1 is None:
        return l2 + '_'
    elif l2 is None:
        return l1 + '_'

    common = []
    for c1, c2 in zip(l1, l2):
        if c1 == c2:
            common.append(c1)
            continue
        break
    common = ''.join(common)
    l1 = l1[len(common):] or '_'
    l2 = l2[len(common):] or '_'
    return '%s-%s,%s' % (common, l1, l2)