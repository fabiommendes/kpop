def test_sharp_classifier(popA, popB):
    pop = popA + popB
    for ind in popA:
        assert pop.cls.parent(ind) == 'A'
    for ind in popB:
        assert pop.cls.parent(ind) == 'B'


def test_probabilistic_classifier(popA, popB):
    pop = popA + popB
    pAs = []

    for ind in popA:
        prob = pop.cls.parent_prob(ind)
        pA = prob['A']
        pAs.append(pA)

    print(pAs)
    assert all([p > 0.75 for p in pAs])