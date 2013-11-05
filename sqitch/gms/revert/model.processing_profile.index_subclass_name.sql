-- Revert model.processing_profile.index_subclass_name

BEGIN;

DROP INDEX model.processing_profile_subclass_index;

COMMIT;
