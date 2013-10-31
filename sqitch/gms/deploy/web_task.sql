-- Deploy web_talk
-- requires: web_schema

BEGIN;

CREATE TABLE IF NOT EXISTS web.task (
    id character varying(255) NOT NULL,
    user_id character varying(255) NOT NULL,
    command_class character varying(255) NOT NULL,
    stdout_pathname character varying(4096),
    stderr_pathname character varying(4096),
    status character varying(50) NOT NULL,
    time_submitted timestamp(6) without time zone DEFAULT now() NOT NULL,
    time_started timestamp(6) without time zone,
    time_finished timestamp(6) without time zone,
    CONSTRAINT task_pkey PRIMARY KEY (id)
);

COMMIT;
