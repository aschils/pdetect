
template<int dim>
PeriodicConstraints<dim>::PeriodicConstraints(){

}

template<int dim>
void PeriodicConstraints<dim>::make_periodicity_constraints(ConstraintMatrix *constraints,
        													DoFHandler<dim> *dof_handler) {
    std::map<unsigned int, double> dof_locations;

    for(DoFHandler<2>::active_cell_iterator cell = dof_handler->begin_active(); 
    			cell != dof_handler->end(); ++cell) {
		if(cell->at_boundary() && cell->face(1)->at_boundary()) {

			dof_locations[cell->face(1)->vertex_dof_index(0, 0)]
            				= cell->face(1)->vertex(0)[1];
          	dof_locations[cell->face(1)->vertex_dof_index(1, 0)]
            				= cell->face(1)->vertex(1)[1];
        }
    }
   
    for (DoFHandler<2>::active_cell_iterator cell = dof_handler->begin_active ();
         			cell != dof_handler->end (); ++cell) {
		if (cell->at_boundary() && cell->face(0)->at_boundary()) {
          	for (unsigned int face_vertex = 0; face_vertex<2; ++face_vertex) {

				constraints->add_line(cell->face(0)->vertex_dof_index(face_vertex, 0));
				std::map<unsigned int, double>::const_iterator p = dof_locations.begin();
              	for(; p != dof_locations.end(); ++p) {
                	if(std::fabs(p->second-cell->face(0)->vertex(face_vertex)[1])
                					< 1e-8) {
						constraints->add_entry(cell->face(0)->vertex_dof_index(face_vertex, 0),
                                    p->first, 1.0);
                    	break;
                  	}
                }
              	Assert (p != dof_locations.end(),
                    ExcMessage ("No corresponding degree of freedom was found!"));
            }
        }
    }
}
