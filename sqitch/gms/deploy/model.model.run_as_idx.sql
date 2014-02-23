-- Deploy model.model.run_as_idx
-- requires: model.model.run_as.sql

CREATE INDEX CONCURRENTLY model__run_as_idx ON model.model (run_as);
