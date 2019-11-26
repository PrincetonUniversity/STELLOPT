@0xcb6e6711bc51954d;
interface Eval {
  eval @0 (params :List(Float64)) -> (objective :Float64);
  evalAll @1 (params :List(Float64)) -> (objectives :List(Float64));
}
