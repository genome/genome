-- Deploy timeline.base.object_id
-- requires: timeline_base

BEGIN;

CREATE INDEX base_object_id_idx on timeline.base using btree (object_id);

COMMIT;
