from .attr import Attr


class Admixture(Attr):
    def structure_admixture(self, k, pop_ids=None, parental_ids=None):
        """
        Runs the ADMIXTURE program to detect structure in the population.

        Args:
            k: number of parental populations.

        Returns:
            A new Population object with all individuals classified with their
            respective admixture coefficients.
        """

        from kpop.external.admixture.base import run_admixture
        pop = self._population

        if not pop.is_biallelic:
            raise ValueError('ADMIXTURE only supports biallelic populations')

        kwargs = {}
        if pop_ids:
            kwargs['supervised'] = pop_ids
        result = run_admixture(pop, k, disp=0, **kwargs)
        parental = result.make_parental(ids=parental_ids)
        individuals = result.make_admixture_ids(pop)
        out = pop.transformed_copy(individuals, parent=parental)
        out.admixture_result = result
        return out
