-- Deploy model.model.created_by_idx
-- requires: model.model.created_by

CREATE INDEX CONCURRENTLY model__created_by_idx ON model.model (created_by);
