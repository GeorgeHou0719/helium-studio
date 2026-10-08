"""Boundary-aware COMSOL profile sampling without requiring a COMSOL license."""
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, MagicMock
import numpy as np
from models.fluid_eda.comsol_worker import PipeSolver


class ProfileSamplingTests(unittest.TestCase):
    def test_exact_endpoints_and_nearby_interior_points_keep_coordinates_and_order(self):
        for length in (.0001, .5, 1.5):
            for region, ends in [('gas', ('inlet', 'outlet')), ('steel', ('steelends', 'steelwarm'))]:
                with self.subTest(length=length, region=region):
                    solver = PipeSolver.__new__(PipeSolver)
                    radii = np.array([.001, .002, .003, .004, .005, .006])
                    z = np.array([length, 0, np.nextafter(length, 0), length / 2,
                                  np.nextafter(0., length), length])
                    original_r, original_z = radii.copy(), z.copy()
                    def sample(tag, expressions, r, x, *, boundary=None):
                        self.assertGreater(len(x), 0)
                        if boundary == ends[0]:
                            self.assertTrue(np.all(x == 0))
                        elif boundary == ends[1]:
                            self.assertTrue(np.all(x == length))
                        else:
                            self.assertIsNone(boundary)
                            self.assertTrue(np.all((x != 0) & (x != length)))
                        # Signed velocity is retained, independently of point order.
                        return np.vstack([r, x, -r-x])
                    solver._sample = Mock(side_effect=sample)
                    actual = solver._sample_profile('profile', ['r','z','w'], radii, z,
                                                    length=length, region=region)
                    np.testing.assert_array_equal(actual, np.vstack([radii,z,-radii-z]))
                    np.testing.assert_array_equal(radii, original_r)
                    np.testing.assert_array_equal(z, original_z)
                    self.assertEqual(solver._sample.call_count, 3)
                    self.assertEqual([c.kwargs.get('boundary') for c in solver._sample.call_args_list],
                                     [None, ends[0], ends[1]])

    def test_single_group_does_not_evaluate_empty_selections(self):
        for z in ([0.,0.], [.5,.5], [.1,.2]):
            solver = PipeSolver.__new__(PipeSolver)
            solver._sample = Mock(return_value=np.array([[3.,4.]]))
            result = solver._sample_profile('one_group', ['T'], [0.,.001], z, length=.5)
            np.testing.assert_array_equal(result, [[3.,4.]])
            solver._sample.assert_called_once()

    def make_solver(self, raw, entities=(27,41)):
        solver = PipeSolver.__new__(PipeSolver)
        j = MagicMock()
        ev = j.result.return_value.numerical.return_value.create.return_value
        ev.getData.return_value = raw
        j.result.return_value.numerical.return_value.tags.return_value = ['test']
        j.component.return_value.selection.return_value.entities.return_value = entities
        solver.model = SimpleNamespace(java=j)
        return solver, j, ev

    def test_boundary_entities_are_resolved_from_current_geometry(self):
        solver, j, ev = self.make_solver(np.array([[[1.,2.]], [[3.,4.]]]))
        actual = solver._sample('test', ['T','w'], [0.,.001], [.5,.5], boundary='outlet')
        np.testing.assert_array_equal(actual, [[1.,2.],[3.,4.]])
        j.component.assert_called_once_with('comp1')
        j.component.return_value.selection.assert_called_once_with('outlet')
        j.component.return_value.selection.return_value.entities.assert_called_once_with(1)
        ev.selection.return_value.geom.assert_called_once_with(1)
        ev.selection.return_value.set.assert_called_once_with([27,41])
        j.result.return_value.numerical.return_value.remove.assert_called_once_with('test')

    def test_invalid_samples_identify_expression_and_location(self):
        solver, _, _ = self.make_solver(np.array([[[1.,2.]], [[3.,np.nan]]]))
        with self.assertRaisesRegex(RuntimeError, 'w at r=0.001 m, z=0.5 m .*outlet'):
            solver._sample('test', ['T','w'], [0.,.001], [.5,.5], boundary='outlet')
        for shape in [(2,2,2), (2,1,1), (2,2)]:
            with self.subTest(shape=shape):
                solver, _, _ = self.make_solver(np.ones(shape))
                with self.assertRaisesRegex(RuntimeError, 'returned shape'):
                    solver._sample('test', ['T','w'], [0.,.001], [.5,.5])

    def test_missing_boundary_does_not_fall_back_to_domain_or_another_entity(self):
        solver, _, ev = self.make_solver(np.ones((1,1,2)), entities=[])
        with self.assertRaisesRegex(RuntimeError, 'boundary selection outlet is empty'):
            solver._sample('test', ['T'], [0.,.001], [.5,.5], boundary='outlet')
        ev.getData.assert_not_called()


if __name__ == '__main__':
    unittest.main()
