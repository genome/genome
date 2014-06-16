-- Revert model.build_status_update_trigger

BEGIN;

DROP TRIGGER build_status_trigger ON model.event;
DROP FUNCTION synchronize_build_status();

COMMIT;
