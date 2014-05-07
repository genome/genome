-- Deploy model.build_status_update_trigger

BEGIN;

CREATE FUNCTION synchronize_build_status() RETURNS trigger AS $synchronize_build_status$
    BEGIN
        UPDATE model.build
        SET status = NEW.event_status, run_by = NEW.user_name, date_scheduled = NEW.date_scheduled, date_completed = NEW.date_completed
        WHERE build_id = NEW.build_id;

        RETURN null;
    END;
$synchronize_build_status$ LANGUAGE plpgsql;

CREATE TRIGGER build_status_trigger
AFTER INSERT OR UPDATE OF event_status,date_scheduled,date_completed,user_name ON model.event
FOR EACH ROW
WHEN (NEW.event_type = 'genome model build')
EXECUTE PROCEDURE synchronize_build_status();

COMMIT;
