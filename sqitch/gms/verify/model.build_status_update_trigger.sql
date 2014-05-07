-- Verify model.build_status_update_trigger

BEGIN;

SELECT 1/COUNT(*) FROM pg_proc WHERE proname = 'synchronize_build_status';

SELECT 1/COUNT(*) FROM pg_trigger WHERE tgname = 'build_status_trigger';

ROLLBACK;
